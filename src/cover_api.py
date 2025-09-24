#!/usr/bin/env python3

"""
FastAPI backend for COVER - COding Variant Effect pRediction
Provides API endpoints for querying transcript databases and finding candidate regions.
"""

import argparse
import logging
import os
import sqlite3
import tempfile
from typing import List, Optional, Union, Dict, Any

import pandas as pd
from fastapi import FastAPI, HTTPException
from pydantic import BaseModel, Field

# Import local modules
from find_candidate_region import main_find_candidate_region
from calculate_het_freq import main_calcualte_het_freq

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='[%(asctime)s] - [%(levelname)s]: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)

# Global database path - will be set when starting the server
DATABASE_PATH: Optional[str] = None
VCF_FILE_PATH: Optional[str] = None

# FastAPI app instance
app = FastAPI(
    title="COVER API",
    description="API for COVER web service",
    version="0.1.0"
)

# Pydantic models for request/response validation
class QueryRequest(BaseModel):
    """Request model for database queries"""
    gene_id: Optional[str] = Field(None, description="Gene ID to search for")
    gene_name: Optional[str] = Field(None, description="Gene name to search for")
    limit: Optional[int] = Field(100, description="Maximum number of results to return")

class TranscriptRecord(BaseModel):
    """Model for transcript record"""
    transcript_id: str
    gene_id: str
    gene_name: str
    seqname: str
    start: int
    end: int
    strand: str
    cds_start: int
    stop_end: int
    exon_list: str

class AnnotationRecord(BaseModel):
    """Model for annotation record"""
    transcript_id: str
    transcript_name: Optional[str]
    transcript_biotype: Optional[str]
    transcript_support_level: Optional[str]
    MANE_select: Optional[bool]
    Canonical: Optional[bool]

class QueryResponse(BaseModel):
    """Response model for database queries"""
    transcripts: List[TranscriptRecord]
    annotations: List[AnnotationRecord]
    total_count: int

class RegionRequest(BaseModel):
    """Request model for candidate region finding"""
    transcript_ids: Union[str, List[str]] = Field(..., description="Single transcript ID or list of transcript IDs")
    max_deletion: Optional[int] = Field(10000, description="Maximum deletion size")
    splice_donor_len: Optional[int] = Field(10, description="Length of splice donor region")
    splice_receptor_len: Optional[int] = Field(28, description="Length of splice receptor region")
    n_before_stop: Optional[int] = Field(2, description="Minimum number of exons before stop codon")

class HetFreqRequest(BaseModel):
    """Request model for heterozygous frequency calculation"""
    regions: List[Dict[str, Any]] = Field(..., description="List of regions in the same format as get_region returns")
    population: str = Field(..., description="Population to analyze (AFR, AMR, EAS, EUR, SAS)")
    max_deletion: Optional[int] = Field(10000, description="Maximum deletion size")
    maf: Optional[float] = Field(0.05, description="MAF cutoff")
    exc_het: Optional[float] = Field(1e-5, description="Excess heterozygosity test p-value cutoff")
    n_pair_max: Optional[int] = Field(200, description="Maximum number of variant pairs to consider")
    pair_het_cutoff: Optional[float] = Field(0.1, description="Minimum heterozygous frequency cutoff")

class HetFreqResponse(BaseModel):
    """Response model for heterozygous frequency results"""
    results: List[Dict[str, Any]]
    total_count: int

class PairHetFreqResponse(BaseModel):
    """Response model for pair heterozygous frequency results"""
    results: List[Dict[str, Any]]
    total_count: int

class PairHetFreqRequest(BaseModel):
    """Request model for pair heterozygous frequency calculation"""
    regions: List[Dict[str, Any]] = Field(..., description="List of regions in the same format as get_region returns")
    population: str = Field(..., description="Population to analyze (AFR, AMR, EAS, EUR, SAS)")
    pair_het_cutoff: float = Field(0.1, description="Minimum heterozygous frequency cutoff for pairs")
    max_deletion: Optional[int] = Field(10000, description="Maximum deletion size")
    maf: Optional[float] = Field(0.05, description="MAF cutoff")
    exc_het: Optional[float] = Field(1e-5, description="Excess heterozygosity test p-value cutoff")
    n_pair_max: Optional[int] = Field(200, description="Maximum number of variant pairs to consider")
    top_n_comb: Optional[int] = Field(1000, description="Top N combinations to output")

class RegionResponse(BaseModel):
    """Response model for candidate regions"""
    results: List[Dict[str, Any]]
    total_count: int
    missing_transcripts: Optional[List[str]] = None

# Dependency to get database connection
def get_db_connection() -> sqlite3.Connection:
    """Get database connection"""
    if DATABASE_PATH is None:
        raise HTTPException(status_code=500, detail="Database not configured")
    
    if not os.path.exists(DATABASE_PATH):
        raise HTTPException(status_code=500, detail=f"Database file not found: {DATABASE_PATH}")
    
    try:
        conn = sqlite3.connect(DATABASE_PATH)
        return conn
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to connect to database: {str(e)}")

# API Endpoints
@app.get("/")
async def root():
    """Root endpoint with API information"""
    return {
        "message": "COVER API",
        "version": "0.1.0",
        "endpoints": {
            "/query_db": "Query database by gene_id or gene_name",
            "/get_region": "Find candidate regions for transcript IDs",
            "/get_het_freq": "Calculate heterozygous frequencies for SNP pairs",
            "/get_pair_het_freq": "Calculate pair heterozygous frequencies for combinations",
            "/health": "Health check endpoint"
        }
    }

@app.get("/health")
async def health_check():
    """Health check endpoint"""
    health_status = {"status": "healthy"}
    issues = []

    try:
        # Check database connectivity
        conn = get_db_connection()
        conn.close()
        health_status["database"] = "connected"
    except Exception as e:
        health_status["database"] = "disconnected"
        issues.append(f"Database error: {str(e)}")

    # Check VCF file configuration and existence
    if VCF_FILE_PATH is None:
        health_status["vcf_file"] = "not configured"
        issues.append("VCF file not configured")
    elif not os.path.exists(VCF_FILE_PATH):
        health_status["vcf_file"] = "not found"
        issues.append(f"VCF file not found: {VCF_FILE_PATH}")
    else:
        health_status["vcf_file"] = "available"

    # Determine overall health status
    if issues:
        health_status["status"] = "unhealthy"
        health_status["issues"] = issues

    return health_status

@app.post("/query_db", response_model=QueryResponse)
async def query_database(request: QueryRequest):
    """
    Query the database for transcripts and annotations by gene_id or gene_name.
    
    Args:
        request: Query parameters including gene_id, gene_name, and limit
    
    Returns:
        QueryResponse containing matching transcripts and annotations
    """
    conn = None
    try:
        # Get database connection for this thread
        conn = get_db_connection()
        
        # Validate input - exactly one search parameter is required
        if not request.gene_id and not request.gene_name:
            raise HTTPException(
                status_code=400, 
                detail="Exactly one of gene_id or gene_name must be provided"
            )
        if request.gene_id and request.gene_name:
            raise HTTPException(
                status_code=400, 
                detail="Cannot specify both gene_id and gene_name. Choose one."
            )
        
        # Build query conditions
        conditions = []
        params = []
        
        if request.gene_id:
            conditions.append("t.gene_id = ?")
            params.append(request.gene_id)
        elif request.gene_name:
            conditions.append("t.gene_name = ?")
            params.append(request.gene_name)
        
        where_clause = " AND ".join(conditions)
        
        # Query transcripts
        transcript_query = f"""
        SELECT t.transcript_id, t.gene_id, t.gene_name, t.seqname, t.start, t.end, 
               t.strand, t.cds_start, t.stop_end, t.exon_list
        FROM transcripts t
        WHERE {where_clause}
        LIMIT ?
        """
        params.append(request.limit)
        
        df_transcripts = pd.read_sql(transcript_query, conn, params=params)
        
        # Query annotations for matching transcript IDs
        if not df_transcripts.empty:
            transcript_ids = df_transcripts['transcript_id'].tolist()
            placeholders = ','.join(['?' for _ in transcript_ids])
            annotation_query = f"""
            SELECT transcript_id, transcript_name, transcript_biotype, 
                   transcript_support_level, MANE_select, Canonical
            FROM annotations
            WHERE transcript_id IN ({placeholders})
            """
            df_annotations = pd.read_sql(annotation_query, conn, params=transcript_ids)
        else:
            df_annotations = pd.DataFrame()
        
        # Convert to response models
        transcripts = [TranscriptRecord(**row) for row in df_transcripts.to_dict('records')]
        annotations = [AnnotationRecord(**row) for row in df_annotations.to_dict('records')]
        
        return QueryResponse(
            transcripts=transcripts,
            annotations=annotations,
            total_count=len(transcripts)
        )
        
    except Exception as e:
        logger.error(f"Error in query_database: {str(e)}")
        raise HTTPException(status_code=500, detail=f"Database query failed: {str(e)}")
    finally:
        # Ensure connection is properly closed in the same thread
        if conn:
            conn.close()

@app.post("/get_region", response_model=RegionResponse)
async def get_candidate_regions(request: RegionRequest):
    """
    Find candidate regions for given transcript ID(s) using main_find_candidate_region.
    
    Args:
        request: Region request parameters including transcript_ids and analysis parameters
    
    Returns:
        RegionResponse containing candidate regions as JSON
    """
    try:
        # Normalize transcript_ids to list
        if isinstance(request.transcript_ids, str):
            transcript_ids = [request.transcript_ids]
        else:
            transcript_ids = request.transcript_ids
        
        if not transcript_ids:
            raise HTTPException(status_code=400, detail="At least one transcript_id must be provided")
        
        # Create temporary files for input and output
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as temp_input:
            temp_input.write('\n'.join(transcript_ids))
            temp_input_path = temp_input.name
        
        with tempfile.TemporaryDirectory() as temp_dir:
            output_prefix = os.path.join(temp_dir, 'candidate_regions')
            
            # Create arguments namespace for main_find_candidate_region
            args = argparse.Namespace(
                id=None,  # Using id_list instead
                id_list=temp_input_path,
                db=DATABASE_PATH,
                output=output_prefix,
                max_deletion=request.max_deletion,
                splice_donor_len=request.splice_donor_len,
                splice_receptor_len=request.splice_receptor_len,
                n_before_stop=request.n_before_stop
            )
            
            # Run the main function
            main_find_candidate_region(args)
            
            # Read the output files
            candidate_file = f"{output_prefix}.candidate_region.txt"
            missing_file = f"{output_prefix}.missing_transcript.txt"
            
            results = []
            missing_transcripts = None
            
            # Read candidate regions if file exists
            if os.path.exists(candidate_file):
                df_results = pd.read_csv(candidate_file, sep='\t')
                results = df_results.to_dict('records')
            
            # Read missing transcripts if file exists
            if os.path.exists(missing_file):
                with open(missing_file, 'r') as f:
                    missing_transcripts = [line.strip() for line in f if line.strip()]
        
        # Clean up temporary input file
        os.unlink(temp_input_path)
        
        return RegionResponse(
            results=results,
            total_count=len(results),
            missing_transcripts=missing_transcripts
        )
        
    except Exception as e:
        logger.error(f"Error in get_candidate_regions: {str(e)}")
        raise HTTPException(status_code=500, detail=f"Region analysis failed: {str(e)}")

@app.post("/get_het_freq", response_model=HetFreqResponse)
async def get_heterozygous_frequencies(request: HetFreqRequest):
    """
    Calculate heterozygous frequencies for SNP pairs in candidate regions.
    This function sets pair_het_cutoff to 1 to skip pair calculations and only return df_het_freq.

    Args:
        request: Request parameters including regions, population, and analysis settings

    Returns:
        HetFreqResponse containing heterozygous frequency results as JSON
    """
    # Validate VCF file is configured
    if VCF_FILE_PATH is None:
        raise HTTPException(status_code=500, detail="VCF file not configured")
    if not os.path.exists(VCF_FILE_PATH):
        raise HTTPException(status_code=500, detail=f"VCF file not found: {VCF_FILE_PATH}")

    try:
        # Create temporary file for regions data
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as temp_regions:
            # Convert regions list to DataFrame and write to temp file
            df_regions = pd.DataFrame(request.regions)
            df_regions.to_csv(temp_regions.name, sep='\t', index=False)
            regions_file = temp_regions.name

        with tempfile.TemporaryDirectory() as temp_dir:
            output_prefix = os.path.join(temp_dir, 'het_freq')

            # Create arguments namespace for main_calculate_het_freq
            args = argparse.Namespace(
                region=regions_file,
                vcf=VCF_FILE_PATH,  # Use global VCF file path
                pop=request.population,
                output=output_prefix,
                index='all',
                maf=request.maf,
                exchet=request.exc_het,
                max_deletion=request.max_deletion,
                n_pair_max=request.n_pair_max,
                pair_het_cutoff=1.0,  # Force to 1 to skip pair calculations
                top_n_comb=1000
            )

            # Run the main function
            main_calcualte_het_freq(args)

            # Read the output file
            het_freq_file = f"{output_prefix}.het_freq.txt"

            results = []
            if os.path.exists(het_freq_file):
                df_results = pd.read_csv(het_freq_file, sep='\t')
                results = df_results.to_dict('records')
            else:
                logger.warning(f"Output file not found: {het_freq_file}")

        # Clean up temporary regions file
        os.unlink(regions_file)

        return HetFreqResponse(
            results=results,
            total_count=len(results)
        )

    except Exception as e:
        logger.error(f"Error in get_heterozygous_frequencies: {str(e)}")
        raise HTTPException(status_code=500, detail=f"Heterozygous frequency calculation failed: {str(e)}")

@app.post("/get_pair_het_freq", response_model=PairHetFreqResponse)
async def get_pair_heterozygous_frequencies(request: PairHetFreqRequest):
    """
    Calculate pair heterozygous frequencies for combinations of SNP pairs in candidate regions.
    This function allows user to specify pair_het_cutoff to control which pairs are considered.

    Args:
        request: Request parameters including regions, population, and analysis settings

    Returns:
        PairHetFreqResponse containing pair heterozygous frequency results as JSON
    """
    # Validate VCF file is configured
    if VCF_FILE_PATH is None:
        raise HTTPException(status_code=500, detail="VCF file not configured")
    if not os.path.exists(VCF_FILE_PATH):
        raise HTTPException(status_code=500, detail=f"VCF file not found: {VCF_FILE_PATH}")

    try:
        # Create temporary file for regions data
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as temp_regions:
            # Convert regions list to DataFrame and write to temp file
            df_regions = pd.DataFrame(request.regions)
            df_regions.to_csv(temp_regions.name, sep='\t', index=False)
            regions_file = temp_regions.name

        with tempfile.TemporaryDirectory() as temp_dir:
            output_prefix = os.path.join(temp_dir, 'pair_het_freq')

            # Create arguments namespace for main_calculate_het_freq
            args = argparse.Namespace(
                region=regions_file,
                vcf=VCF_FILE_PATH,  # Use global VCF file path
                pop=request.population,
                output=output_prefix,
                index='all',
                maf=request.maf,
                exchet=request.exc_het,
                max_deletion=request.max_deletion,
                n_pair_max=request.n_pair_max,
                pair_het_cutoff=request.pair_het_cutoff,  # Use user-specified cutoff
                top_n_comb=request.top_n_comb
            )

            # Run the main function
            main_calcualte_het_freq(args)

            # Read the output file (pair heterozygous frequencies)
            pair_het_freq_file = f"{output_prefix}.pair_het_freq.txt"

            results = []
            if os.path.exists(pair_het_freq_file):
                df_results = pd.read_csv(pair_het_freq_file, sep='\t')
                results = df_results.to_dict('records')
            else:
                logger.warning(f"Output file not found: {pair_het_freq_file}")

        # Clean up temporary regions file
        os.unlink(regions_file)

        return PairHetFreqResponse(
            results=results,
            total_count=len(results)
        )

    except Exception as e:
        logger.error(f"Error in get_pair_heterozygous_frequencies: {str(e)}")
        raise HTTPException(status_code=500, detail=f"Pair heterozygous frequency calculation failed: {str(e)}")

# Function to set database path (called when starting server)
def set_database_path(db_path: str):
    """Set the global database path"""
    global DATABASE_PATH
    DATABASE_PATH = db_path
    logger.info(f"Database path set to: {db_path}")

# Function to set VCF file path (called when starting server)
def set_vcf_file_path(vcf_path: str):
    """Set the global VCF file path"""
    global VCF_FILE_PATH
    VCF_FILE_PATH = vcf_path
    logger.info(f"VCF file path set to: {vcf_path}")

# CLI interface for starting the server
def main():
    """Main function to start the FastAPI server"""
    import uvicorn
    
    parser = argparse.ArgumentParser(description="Start COVER API server")
    parser.add_argument("-d", "--database", required=True, help="Path to SQLite database file")
    parser.add_argument("-v", "--vcf", required=True, help="Path to VCF file")
    parser.add_argument("-p", "--port", type=int, default=8000, help="Port to run the server on")
    parser.add_argument("--host", default="127.0.0.1", help="Host to run the server on")
    
    args = parser.parse_args()
    
    # Validate database file exists
    if not os.path.exists(args.database):
        logger.error(f"Database file not found: {args.database}")
        exit(1)
    
    # Validate VCF file exists
    if not os.path.exists(args.vcf):
        logger.error(f"VCF file not found: {args.vcf}")
        exit(1)
    
    # Set database and VCF paths
    set_database_path(args.database)
    set_vcf_file_path(args.vcf)
    
    # Start server
    logger.info(f"Starting COVER API server on {args.host}:{args.port}")
    logger.info(f"Using database: {args.database}")
    logger.info(f"Using VCF file: {args.vcf}")
    
    uvicorn.run(app, host=args.host, port=args.port)

if __name__ == "__main__":
    main()