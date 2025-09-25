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

# Global database paths - will be set when starting the server
TX_DATABASE_PATH: Optional[str] = None  # Transcript database
RESULTS_DATABASE_PATH: Optional[str] = None  # Results database with het_freq table
VCF_FILE_PATH: Optional[str] = None

# FastAPI app instance
app = FastAPI(
    title="COVER API",
    description="API for COVER web service",
    version="0.1.0"
)

class HetFreqRecord(BaseModel):
    """Model for heterozygous frequency record"""
    transcript_id: str
    gene_id: str
    gene_name: str
    variant1: str
    variant1_region: str
    variant2: str
    variant2_region: str
    distance: int
    target: str
    consequence: str
    cis_het_freq: float
    trans_het_freq: float
    max_het_freq: float
    target_genotype: str

class ResultsRequest(BaseModel):
    """Request model for querying analysis results"""
    gene_id: Optional[str] = Field(None, description="Gene ID to filter by")
    transcript_id: Optional[str] = Field(None, description="Transcript ID to filter by")
    variant_id: Optional[str] = Field(None, description="Variant ID to filter by (searches both variant1 and variant2)")
    limit: Optional[int] = Field(100, description="Maximum number of results to return")

class ResultsResponse(BaseModel):
    """Response model for analysis results queries"""
    results: List[HetFreqRecord]
    total_count: int

# Pydantic models for request/response validation
class TranscriptRequest(BaseModel):
    """Request model for transcript queries"""
    gene_id: Optional[str] = Field(None, description="Gene ID to search for")
    gene_name: Optional[str] = Field(None, description="Gene name to search for")
    limit: Optional[int] = Field(100, description="Maximum number of results to return")

class TranscriptRecord(BaseModel):
    """Model for transcript record"""
    transcript_id: str
    gene_id: str
    gene_name: str
    transcript_name: Optional[str]
    transcript_biotype: Optional[str]
    transcript_support_level: Optional[str]
    MANE_select: Optional[bool]
    Canonical: Optional[bool]

class ExonRecord(BaseModel):
    """Model for exon record"""
    exon_id: str
    transcript_id: str
    seqname: str
    start: int
    end: int
    strand: str

class TranscriptResponse(BaseModel):
    """Response model for transcript queries"""
    transcript_table: List[TranscriptRecord]
    exon_table: List[ExonRecord]
    total_transcripts: int
    total_exons: int

class RegionRequest(BaseModel):
    """Request model for candidate region finding"""
    transcript_ids: Union[str, List[str]] = Field(..., description="Single transcript ID or list of transcript IDs")
    max_deletion: Optional[int] = Field(10000, description="Maximum deletion size")
    splice_donor_len: Optional[int] = Field(10, description="Length of splice donor region")
    splice_receptor_len: Optional[int] = Field(28, description="Length of splice receptor region")
    n_before_stop: Optional[int] = Field(2, description="Minimum number of exons before stop codon")

class RegionRecord(BaseModel):
    """Model for region record"""
    transcript_id: str
    gene_id: str
    gene_name: str
    seqname: str
    upstream: str
    upstream_start: int
    upstream_end: int
    downstream: str
    downstream_start: int
    downstream_end: int
    distance: int
    strand: str
    target_exon: str
    consequence: str

class PairHetFreqRecord(BaseModel):
    """Model for pair heterozygous frequency record"""
    transcript_id: str
    gene_id: str
    gene_name: str
    pair1: str
    target1: str
    pair2: str
    target2: str
    pair_het_freq: float

class HetFreqRequest(BaseModel):
    """Request model for heterozygous frequency calculation"""
    regions: List[Dict[str, Any]] = Field(..., description="List of regions in the same format as get_region returns")
    population: str = Field(..., description="Population to analyze (AFR, AMR, EAS, EUR, SAS)")
    max_deletion: Optional[int] = Field(10000, description="Maximum deletion size")
    maf: Optional[float] = Field(0.05, description="MAF cutoff")
    exc_het: Optional[float] = Field(1e-5, description="Excess heterozygosity test p-value cutoff")
    n_pair_max: Optional[int] = Field(200, description="Maximum number of variant pairs to consider")

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
    results: List[RegionRecord]
    total_count: int
    missing_transcripts: Optional[List[str]] = None

class HetFreqResponse(BaseModel):
    """Response model for heterozygous frequency results"""
    results: List[HetFreqRecord]
    total_count: int
    population: str

class PairHetFreqResponse(BaseModel):
    """Response model for pair heterozygous frequency results"""
    results: List[PairHetFreqRecord]
    total_count: int

# Dependency to get database connection
def get_db_connection() -> sqlite3.Connection:
    """Get database connection"""
    if TX_DATABASE_PATH is None:
        raise HTTPException(status_code=500, detail="Database not configured")
    
    if not os.path.exists(TX_DATABASE_PATH):
        raise HTTPException(status_code=500, detail=f"Database file not found: {TX_DATABASE_PATH}")
    
    try:
        conn = sqlite3.connect(TX_DATABASE_PATH)
        return conn
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to connect to database: {str(e)}")

# Dependency to get results database connection
def get_results_db_connection() -> sqlite3.Connection:
    """Get results database connection"""
    if RESULTS_DATABASE_PATH is None:
        raise HTTPException(status_code=500, detail="Results database not configured")

    if not os.path.exists(RESULTS_DATABASE_PATH):
        raise HTTPException(status_code=500, detail=f"Results database file not found: {RESULTS_DATABASE_PATH}")

    try:
        conn = sqlite3.connect(RESULTS_DATABASE_PATH)
        return conn
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to connect to results database: {str(e)}")

# API Endpoints
@app.get("/")
async def root():
    """Root endpoint with API information"""
    return {
        "message": "COVER API",
        "version": "0.1.0",
        "endpoints": {
            "/query_transcript": "Query database by gene_id or gene_name (TranscriptRequest/TranscriptResponse)",
            "/query_results": "Query het_freq results by gene_id, transcript_id, or variant_id (ResultsRequest/ResultsResponse)",
            "/get_region": "Find candidate regions for transcript IDs (RegionRequest/RegionResponse)",
            "/get_het_freq": "Calculate heterozygous frequencies for SNP pairs (HetFreqRequest/HetFreqResponse)",
            "/get_pair_het_freq": "Calculate pair heterozygous frequencies for combinations (PairHetFreqRequest/PairHetFreqResponse)",
            "/health": "Health check endpoint"
        }
    }

@app.get("/health")
async def health_check():
    """Health check endpoint"""
    health_status = {"status": "healthy"}
    issues = []

    try:
        # Check main database connectivity
        conn = get_db_connection()
        conn.close()
        health_status["database"] = "connected"
    except Exception as e:
        health_status["database"] = "disconnected"
        issues.append(f"Database error: {str(e)}")

    # Check results database connectivity
    try:
        conn = get_results_db_connection()
        conn.close()
        health_status["results_database"] = "connected"
    except Exception as e:
        health_status["results_database"] = "disconnected"
        issues.append(f"Results database error: {str(e)}")

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

@app.post("/query_results", response_model=ResultsResponse)
async def query_results(request: ResultsRequest):
    """
    Query the results database for heterozygous frequency data.

    Supports filtering by gene_id, transcript_id, or variant_id (searches both variant1 and variant2).
    At least one filter parameter must be provided.

    Example usage:
        "gene_id": "ENSG00000142192", "transcript_id": "ENST00000354192", "variant_id": "chr21:26161857:T:G"

    Returns:
        ResultsResponse containing filtered heterozygous frequency results
    """
    conn = None
    try:
        # Get results database connection
        conn = get_results_db_connection()

        # Validate input - at least one search parameter is required
        if not request.gene_id and not request.transcript_id and not request.variant_id:
            raise HTTPException(
                status_code=400,
                detail="At least one of gene_id, transcript_id, or variant_id must be provided"
            )

        # Build query conditions
        conditions = []
        params = []

        if request.gene_id:
            conditions.append("gene_id = ?")
            params.append(request.gene_id)

        if request.transcript_id:
            conditions.append("transcript_id = ?")
            params.append(request.transcript_id)

        if request.variant_id:
            conditions.append("(variant1 = ? OR variant2 = ?)")
            params.extend([request.variant_id, request.variant_id])

        where_clause = " AND ".join(conditions)

        # Query het_freq data with filters
        query = f"""
        SELECT transcript_id, gene_id, gene_name, variant1, variant1_region,
               variant2, variant2_region, distance, target, consequence,
               cis_het_freq, trans_het_freq, max_het_freq, target_genotype
        FROM het_freq
        WHERE {where_clause}
        LIMIT ?
        """
        params.append(request.limit)

        df_results = pd.read_sql(query, conn, params=params)

        # Convert to response models - use vectorized operations for better performance
        if not df_results.empty:
            # Convert data types in bulk using vectorized operations
            df_results = df_results.astype({
                'transcript_id': 'string',
                'gene_id': 'string',
                'gene_name': 'string',
                'variant1': 'string',
                'variant1_region': 'string',
                'variant2': 'string',
                'variant2_region': 'string',
                'distance': 'int64',
                'target': 'string',
                'consequence': 'string',
                'cis_het_freq': 'float64',
                'trans_het_freq': 'float64',
                'max_het_freq': 'float64',
                'target_genotype': 'string'
            })

            # Convert to list of dictionaries using pandas' optimized method
            results = [HetFreqRecord(**row) for row in df_results.to_dict('records')]
        else:
            results = []

        return ResultsResponse(
            results=results,
            total_count=len(results)
        )

    except Exception as e:
        logger.error(f"Error in query_het_freq: {str(e)}")
        raise HTTPException(status_code=500, detail=f"HetFreq query failed: {str(e)}")
    finally:
        # Ensure connection is properly closed
        if conn:
            conn.close()

@app.post("/step1", response_model=TranscriptResponse)
async def query_transcript(request: TranscriptRequest):
    """
    Query the database for transcripts and annotations by gene_id or gene_name.

    Example usage:
        "gene_id": "ENSG00000142192", "gene_name": "APP"

    Returns:
        TranscriptResponse containing merged transcript table and exon table
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

        # Query merged transcript and annotation data (including exon_list for processing)
        transcript_query = f"""
        SELECT t.transcript_id, t.gene_id, t.gene_name, t.exon_list,
               a.transcript_name, a.transcript_biotype, a.transcript_support_level,
               a.MANE_select, a.Canonical
        FROM transcripts t
        LEFT JOIN annotations a ON t.transcript_id = a.transcript_id
        WHERE {where_clause}
        LIMIT ?
        """
        params.append(request.limit)

        df_transcripts = pd.read_sql(transcript_query, conn, params=params)

        # Query all exons for matched transcripts with transcript_id mapping
        if not df_transcripts.empty:
            # Extract all unique exon_ids from exon_list column
            exon_ids = []
            transcript_to_exons = {}

            for _, row in df_transcripts.iterrows():
                transcript_id = row['transcript_id']
                if pd.notna(row['exon_list']) and row['exon_list']:
                    exons = row['exon_list'].split(',')
                    exon_ids.extend(exons)
                    transcript_to_exons[transcript_id] = exons

            if exon_ids:
                exon_placeholders = ','.join(['?' for _ in exon_ids])
                exon_query = f"""
                SELECT exon_id, seqname, start, end, strand
                FROM exons
                WHERE exon_id IN ({exon_placeholders})
                """
                df_exons = pd.read_sql(exon_query, conn, params=exon_ids)

                # Add transcript_id to each exon by mapping back
                df_exons['transcript_id'] = None
                for transcript_id, exons in transcript_to_exons.items():
                    for exon_id in exons:
                        df_exons.loc[df_exons['exon_id'] == exon_id, 'transcript_id'] = transcript_id
            else:
                df_exons = pd.DataFrame()
        else:
            df_exons = pd.DataFrame()

        # Remove exon_list from transcript data before returning
        if not df_transcripts.empty:
            df_transcripts = df_transcripts.drop(columns=['exon_list'])

        # Convert to response models
        transcript_table = [TranscriptRecord(**row) for row in df_transcripts.to_dict('records')]
        exon_table = [ExonRecord(**row) for row in df_exons.to_dict('records')]

        return TranscriptResponse(
            transcript_table=transcript_table,
            exon_table=exon_table,
            total_transcripts=len(transcript_table),
            total_exons=len(exon_table)
        )

    except Exception as e:
        logger.error(f"Error in query_database: {str(e)}")
        raise HTTPException(status_code=500, detail=f"Database query failed: {str(e)}")
    finally:
        # Ensure connection is properly closed in the same thread
        if conn:
            conn.close()

@app.post("/step2", response_model=RegionResponse)
async def get_candidate_regions(request: RegionRequest):
    """
    Find candidate regions for given transcript_id(s).
    
    Example usage:
        "transcript_ids": ["ENST00000354192", "ENST00000348990"]
    
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
                db=TX_DATABASE_PATH,
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
                results = [RegionRecord(**row) for row in df_results.to_dict('records')]
            
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

@app.post("/step3", response_model=HetFreqResponse)
async def get_heterozygous_frequencies(request: HetFreqRequest):
    """
    Calculate heterozygous frequencies for SNP pairs in candidate regions.
    This function sets pair_het_cutoff to 1 to skip pair calculations and only return df_het_freq.

    Example usage:
        "regions": [
            {
            "transcript_id": "ENST00000354192",
            "gene_id": "ENSG00000142192",
            "gene_name": "APP",
            "seqname": "chr21",
            "upstream": "5' region",
            "upstream_start": 26170619,
            "upstream_end": 26180618,
            "downstream": "intron 1",
            "downstream_start": 26090101,
            "downstream_end": 26170553,
            "distance": 66,
            "strand": "-",
            "target_exon": "exon 1",
            "consequence": "start loss"
            }
        ],
        "population": "EUR"


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
                # Add transcript_id, gene_id, and gene_name from the first region in the request
                first_region = request.regions[0] if request.regions else {}
                transcript_id = first_region.get('transcript_id', '')
                gene_id = first_region.get('gene_id', '')
                gene_name = first_region.get('gene_name', '')
                
                # Add these fields to each row before creating records
                df_results['transcript_id'] = transcript_id
                df_results['gene_id'] = gene_id
                df_results['gene_name'] = gene_name
                
                # Remove population from records since it's now in the response
                if 'population' in df_results.columns:
                    df_results = df_results.drop(columns=['population'])
                
                results = [HetFreqRecord(**row) for row in df_results.to_dict('records')]
            else:
                logger.warning(f"Output file not found: {het_freq_file}")

        # Clean up temporary regions file
        os.unlink(regions_file)

        return HetFreqResponse(
            results=results,
            total_count=len(results),
            population=request.population
        )

    except Exception as e:
        logger.error(f"Error in get_heterozygous_frequencies: {str(e)}")
        raise HTTPException(status_code=500, detail=f"Heterozygous frequency calculation failed: {str(e)}")

@app.post("/step4", response_model=PairHetFreqResponse)
async def get_pair_heterozygous_frequencies(request: PairHetFreqRequest):
    """
    Calculate pair heterozygous frequencies for combinations of SNP pairs in candidate regions.
    This function allows user to specify pair_het_cutoff to control which pairs are considered.

    Example usage:
        "regions": [
            {
            "transcript_id": "ENST00000354192",
            "gene_id": "ENSG00000142192",
            "gene_name": "APP",
            "seqname": "chr21",
            "upstream": "5' region",
            "upstream_start": 26170619,
            "upstream_end": 26180618,
            "downstream": "intron 1",
            "downstream_start": 26090101,
            "downstream_end": 26170553,
            "distance": 66,
            "strand": "-",
            "target_exon": "exon 1",
            "consequence": "start loss"
            }
        ],
        "population": "EUR"

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
                # Add transcript_id, gene_id, and gene_name from the first region in the request
                first_region = request.regions[0] if request.regions else {}
                transcript_id = first_region.get('transcript_id', '')
                gene_id = first_region.get('gene_id', '')
                gene_name = first_region.get('gene_name', '')
                
                # Add these fields to each row before creating records
                df_results['transcript_id'] = transcript_id
                df_results['gene_id'] = gene_id
                df_results['gene_name'] = gene_name
                
                results = [PairHetFreqRecord(**row) for row in df_results.to_dict('records')]
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

# Function to set results database path (called when starting server)
def set_results_database_path(db_path: str):
    """Set the global results database path"""
    global RESULTS_DATABASE_PATH
    RESULTS_DATABASE_PATH = db_path
    logger.info(f"Results database path set to: {db_path}")

# Function to set transcript database path (called when starting server)
def set_transcript_database_path(db_path: str):
    """Set the global transcript database path"""
    global TX_DATABASE_PATH
    TX_DATABASE_PATH = db_path
    logger.info(f"Transcript database path set to: {db_path}")

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
    parser.add_argument("-t", "--transcript-database", required=True, help="Path to transcript SQLite database file")
    parser.add_argument("-r", "--results-database", required=True, help="Path to results SQLite database file (with het_freq table)")
    parser.add_argument("-v", "--vcf", required=True, help="Path to VCF file")
    parser.add_argument("-p", "--port", type=int, default=8000, help="Port to run the server on")
    parser.add_argument("--host", default="127.0.0.1", help="Host to run the server on")

    args = parser.parse_args()

    # Validate database file exists
    if not os.path.exists(args.transcript_database):
        logger.error(f"Transcript database file not found: {args.transcript_database}")
        exit(1)

    # Validate results database file exists
    if not os.path.exists(args.results_database):
        logger.error(f"Results database file not found: {args.results_database}")
        exit(1)

    # Validate VCF file exists
    if not os.path.exists(args.vcf):
        logger.error(f"VCF file not found: {args.vcf}")
        exit(1)

    # Set database and VCF paths
    set_transcript_database_path(args.transcript_database)
    set_results_database_path(args.results_database)
    set_vcf_file_path(args.vcf)

    # Start server
    logger.info(f"Starting COVER API server on {args.host}:{args.port}")
    logger.info(f"Using transcript database: {args.transcript_database}")
    logger.info(f"Using results database: {args.results_database}")
    logger.info(f"Using VCF file: {args.vcf}")

    uvicorn.run(app, host=args.host, port=args.port)

if __name__ == "__main__":
    main()