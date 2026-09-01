#!/usr/bin/env python3

"""
FastAPI backend for COVER.
"""

import argparse
import hashlib
import json
import logging
import os
import shutil
import sqlite3
import tempfile
import time
import datetime
from typing import List, Optional, Union, Dict, Any, Tuple

import pandas as pd
from fastapi import FastAPI, HTTPException
from pydantic import BaseModel, Field

# Import local modules
from find_candidate_region import main_find_candidate_region
from calculate_het_freq import main_calculate_het_freq

# Set up logging to both console and file
def setup_logging():
    """Set up logging configuration for both console and file output"""
    # Create logs directory if it doesn't exist
    logs_dir = os.path.join(os.getcwd(), 'logs')
    os.makedirs(logs_dir, exist_ok=True)

    # Create log filename with current date
    current_date = datetime.datetime.now().strftime('%Y-%m-%d')
    log_filename = f'cover-api.{current_date}.log'
    log_filepath = os.path.join(logs_dir, log_filename)

    # Create logger
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)

    # Remove existing handlers to avoid duplicates
    for handler in logger.handlers[:]:
        logger.removeHandler(handler)

    # Create formatters
    file_formatter = logging.Formatter(
        '[%(asctime)s] - [%(levelname)s] - [%(name)s:%(lineno)d]: %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    console_formatter = logging.Formatter(
        '[%(asctime)s] - [%(levelname)s]: %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )

    # File handler for detailed logging
    file_handler = logging.FileHandler(log_filepath, encoding='utf-8')
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(file_formatter)

    # Console handler for basic logging
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(console_formatter)

    # Add handlers to logger
    logger.addHandler(file_handler)
    logger.addHandler(console_handler)

    return logger

# Initialize logging with default level
logger = setup_logging()

# Global database and VCF file paths
TX_DATABASE_PATH = None
RESULTS_DATABASE_PATH = None
VCF_FILE_PATH = None
# Cache directory for storing temporary results
CACHE_DIR = os.path.join(tempfile.gettempdir(), 'cover_cache')
os.makedirs(CACHE_DIR, exist_ok=True)

# Cache expiration time in seconds (1 week)
CACHE_EXPIRATION = 604800

# Global counter for cache cleanup
CACHE_CLEANUP_COUNTER = 0
CLEANUP_THRESHOLD = 1000

def create_cache_key(params: Dict[str, Any], endpoint_name: str) -> str:
    """Create a cache key from parameters dictionary using hash for filesystem safety"""
    # Sort parameters to ensure consistent cache keys
    sorted_params = sorted(params.items())

    # Create a string representation
    param_str = f"{endpoint_name}:{','.join(f'{k}={v}' for k, v in sorted_params)}"

    # Create hash for filesystem-safe key
    cache_hash = hashlib.md5(param_str.encode()).hexdigest()
    return cache_hash

def get_cache_path(cache_key: str, endpoint_name: str) -> str:
    """Get the cache directory path for a given cache key and endpoint"""
    return os.path.join(CACHE_DIR, endpoint_name, cache_key)

def is_cache_valid(cache_path: str, result_files: List[str] = None) -> bool:
    """Check if cache exists, is not expired, and contains required result files"""
    if not os.path.exists(cache_path):
        return False

    # Check if cache is expired
    cache_age = time.time() - os.path.getmtime(cache_path)
    if cache_age >= CACHE_EXPIRATION:
        return False

    # If specific result files are specified, check that they exist
    if result_files:
        for filename in result_files:
            result_file_path = os.path.join(cache_path, filename)
            if not os.path.exists(result_file_path):
                return False

    return True

def load_cached_results(cache_path: str, result_files: Dict[str, str]) -> Dict[str, Any]:
    """Load cached results from files"""
    results = {}

    for file_key, filename in result_files.items():
        file_path = os.path.join(cache_path, filename)
        logger.debug(f"Loading cached file: {file_path} -> exists: {os.path.exists(file_path)}")
        if os.path.exists(file_path):
            # Handle missing_transcript.txt files differently (text format, no header)
            if filename.endswith('missing_transcript.txt'):
                try:
                    with open(file_path, 'r') as f:
                        content = f.read().strip()
                        if content:
                            results[file_key] = [line.strip() for line in content.split('\n') if line.strip()]
                        else:
                            results[file_key] = []  # Empty list for empty file
                    logger.debug(f"Loaded {len(results[file_key])} missing transcripts from {filename}")
                except Exception as e:
                    logger.warning(f"Error reading missing transcript file {file_path}: {e}")
                    results[file_key] = []
            else:
                # All other cached files are TSV format, read with pandas
                try:
                    df = pd.read_csv(file_path, sep='\t')
                    results[file_key] = df
                    logger.debug(f"Loaded {len(df)} records from {filename}")
                except Exception as e:
                    logger.warning(f"Error parsing cached CSV file {file_path}: {e}")
                    results[file_key] = pd.DataFrame()
        else:
            # All cached files should exist, so warn if any are missing
            logger.warning(f"Cached file not found: {file_path}")
            results[file_key] = None

    return results

def initialize_pagination_params(page_limit: Optional[int] = None, page_no: Optional[int] = None) -> Tuple[int, int]:
    """
    Initialize pagination parameters
    """
    
    # Set default values for pagination
    page_limit = page_limit if page_limit is not None else 20
    page_no = page_no if page_no is not None else 1

    # Ensure page_limit is positive
    if page_limit <= 0:
        raise HTTPException(status_code=400, detail="page_limit must be greater than 0")

    # Ensure page_no is positive
    if page_no <= 0:
        raise HTTPException(status_code=400, detail="page_no must be greater than 0")
    
    return page_limit, page_no

def get_pagination_info(all_results: Any, page_limit: Optional[int] = None, page_no: Optional[int] = None) -> Dict[str, Any]:
    """
    Calculate pagination information and return paginated results.

    Args:
        all_results: List of all results
        page_limit: Maximum number of results per page (default: 20)
        page_no: Page number to retrieve (1-indexed, default: 1)

    Returns:
        Dict containing pagination info and paginated results

    Raises:
        HTTPException: If pagination parameters are invalid
    """

    # Support pandas DataFrame and list-like inputs
    if isinstance(all_results, pd.DataFrame):
        total_count = len(all_results)
        total_pages = (total_count + page_limit - 1) // page_limit  # Ceiling division
        start_idx = (page_no - 1) * page_limit
        end_idx = start_idx + page_limit
        paginated_df = all_results.iloc[start_idx:end_idx]
        paginated_results = paginated_df.to_dict('records')
        return {
            'results': paginated_results,
            'total_count': total_count,
            'page_no': page_no,
            'page_limit': page_limit,
            'total_pages': total_pages
        }
    else:
        total_count = len(all_results) if all_results is not None else 0
        total_pages = (total_count + page_limit - 1) // page_limit  # Ceiling division
        start_idx = (page_no - 1) * page_limit
        end_idx = start_idx + page_limit
        if total_count > 0:
            paginated_results = all_results[start_idx:end_idx]
        else:
            paginated_results = []
        return {
            'results': paginated_results,
            'total_count': total_count,
            'page_no': page_no,
            'page_limit': page_limit,
            'total_pages': total_pages
        }

def save_input_parameters(request: Any, cache_path: str) -> None:
    """Save the full input parameters to a JSON file in the cache directory"""
    try:
        # Ensure cache directory exists
        os.makedirs(cache_path, exist_ok=True)

        # Create parameters.json file path
        params_file = os.path.join(cache_path, 'parameters.json')

        # Skip saving if parameters file already exists
        if os.path.exists(params_file):
            logger.debug(f"Parameters file already exists, skipping save: {params_file}")
            return

        # Convert request object to dictionary, excluding pagination and filter parameters
        request_dict = request.model_dump()

        # Remove pagination parameters as they don't affect computation results
        request_dict.pop('page_limit', None)
        request_dict.pop('page_no', None)

        # Remove filter parameters as they are applied after computation/caching
        request_dict.pop('filter_min', None)
        request_dict.pop('filter_max', None)

        # Save to JSON file with indentation for readability
        with open(params_file, 'w', encoding='utf-8') as f:
            json.dump(request_dict, f, indent=2, ensure_ascii=False)

        logger.debug(f"Saved input parameters to: {params_file}")
    except Exception as e:
        logger.warning(f"Failed to save input parameters to {cache_path}: {str(e)}")

def cleanup_expired_cache():
    """Clean up expired cache directories when counter reaches threshold"""
    global CACHE_CLEANUP_COUNTER
    CACHE_CLEANUP_COUNTER = 0  # Reset counter after cleanup

    try:
        current_time = time.time()
        for endpoint_dir in os.listdir(CACHE_DIR):
            endpoint_path = os.path.join(CACHE_DIR, endpoint_dir)
            if os.path.isdir(endpoint_path):
                for cache_item in os.listdir(endpoint_path):
                    cache_item_path = os.path.join(endpoint_path, cache_item)
                    if os.path.isdir(cache_item_path):
                        item_age = current_time - os.path.getmtime(cache_item_path)
                        if item_age > CACHE_EXPIRATION:
                            shutil.rmtree(cache_item_path)
                            logger.info(f"Cleaned up expired cache: {endpoint_dir}/{cache_item}")
    except Exception as e:
        logger.warning(f"Error cleaning up cache: {str(e)}")

def increment_cleanup_counter():
    """Increment cleanup counter and trigger cleanup if threshold reached"""
    global CACHE_CLEANUP_COUNTER
    CACHE_CLEANUP_COUNTER += 1

    if CACHE_CLEANUP_COUNTER >= CLEANUP_THRESHOLD:
        cleanup_expired_cache()

def add_het_freq_meta(df: pd.DataFrame, regions: List[Dict[str, Any]], population: str) -> pd.DataFrame:
    """Add region information to the result DataFrame"""

    # Add transcript_id, gene_id, gene_name, and population from the first region in the request
    first_region = regions[0] if regions else {}
    transcript_id = first_region.get('transcript_id', '')
    gene_id = first_region.get('gene_id', '')
    gene_name = first_region.get('gene_name', '')

    # Add these fields to each row
    df['transcript_id'] = transcript_id
    df['gene_id'] = gene_id
    df['gene_name'] = gene_name
    df['population'] = population

    return df

def run_het_freq_calculation(
    regions: List[Dict[str, Any]],
    population: str,
    cache_key: str,
    cache_path: str,
    pair_het_cutoff: float,
    top_n_comb: int = 1000,
    output_file_suffix: str = "het_freq.het_freq.txt",
    process_pair_results: bool = False
) -> pd.DataFrame:
    """
    Shared function to process regions and calculate heterozygous frequencies.

    Args:
        regions: List of region dictionaries
        population: Population to analyze
        cache_key: Cache key for this operation
        cache_path: Cache path for this operation
        pair_het_cutoff: Cutoff for pair heterozygous frequency
        top_n_comb: Top N combinations to output (default: 1000)
        output_file_suffix: Suffix for output file to read
        process_pair_results: Whether to process pair results (step4) or single results (step3)

    Returns:
        pandas DataFrame of processed results
    """
    logger.info(f"Computing heterozygous frequencies for cache key: {cache_key}")

    # Create temporary file for regions data
    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as temp_regions:
        # Convert regions list to DataFrame and write to temp file
        df_regions = pd.DataFrame(regions)
        df_regions.to_csv(temp_regions.name, sep='\t', index=False)
        regions_file = temp_regions.name

    try:
        # Write directly to cache directory
        output_prefix = os.path.join(cache_path, 'het_freq' if not process_pair_results else 'pair_het_freq')

        # Create arguments namespace for main_calculate_het_freq
        args = argparse.Namespace(
            region=regions_file,
            vcf=VCF_FILE_PATH,  # Use global VCF file path
            pop=population,
            output=output_prefix,
            index='all',
            maf=0.05,  # Default MAF
            exchet=1e-5,  # Default excess heterozygosity cutoff
            max_deletion=10000,  # Default max deletion
            n_pair_max=200,  # Default n_pair_max
            pair_het_cutoff=pair_het_cutoff,
            top_n_comb=top_n_comb
        )

        # Run the main function
        main_calculate_het_freq(args)

        # Read the output file
        output_file = f"{output_prefix}.{output_file_suffix}"

        df_results = pd.DataFrame()
        if os.path.exists(output_file):
            # if no suitable snps, the file is empty
            try:
                df_results = pd.read_csv(output_file, sep='\t')
                df_results = add_het_freq_meta(df_results, regions, population)
            except pd.errors.EmptyDataError:
                logger.warning(f"Empty step3 output for cache key: {cache_key}")
                return df_results
        else:
            logger.error(f"Output file not found: {output_file}")

    finally:
        # Clean up temporary regions file
        os.unlink(regions_file)

    return df_results

# FastAPI app instance
app = FastAPI(
    title="COVER API",
    description="API for COVER web service",
    version="0.1.0"
)

# Configure CORS
def configure_cors(app: FastAPI, allowed_origins: List[str]) -> None:
    from fastapi.middleware.cors import CORSMiddleware
    app.add_middleware(
        CORSMiddleware,
        allow_origins=allowed_origins,
        allow_credentials=True,
        allow_methods=["*"],  # Allows all methods
        allow_headers=["*"],  # Allows all headers
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
    population: str

class ResultsRequest(BaseModel):
    """Request model for querying analysis results"""
    gene_id: Optional[str] = Field("ENSG00000142192", description="Gene ID to filter by")
    transcript_id: Optional[str] = Field(None, description="Transcript ID to filter by")
    variant_id: Optional[str] = Field(None, description="Variant ID to filter by (searches both variant1 and variant2)")
    page_limit: Optional[int] = Field(20, description="Maximum number of results per page")
    page_no: Optional[int] = Field(1, description="Page number to retrieve (starting from 1)")

class ResultsResponse(BaseModel):
    """Response model for analysis results queries"""
    results: List[HetFreqRecord]
    total_count: int
    page_no: int
    page_limit: int
    total_pages: int

# Pydantic models for request/response validation
class TranscriptRequest(BaseModel):
    """Request model for transcript queries"""
    gene_id: Optional[str] = Field(None, description="Gene ID to search for")
    gene_name: Optional[str] = Field(None, description="Gene name to search for")
    limit: Optional[int] = Field(100, description="Maximum number of results to return")

class TranscriptRecord(BaseModel):
    """Model for transcript record"""
    transcript_no: int
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
    exon_no: int
    transcript_id: str
    start: int
    end: int

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
    page_limit: Optional[int] = Field(20, description="Maximum number of results per page")
    page_no: Optional[int] = Field(1, description="Page number to retrieve (starting from 1)")

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
    page_limit: Optional[int] = Field(20, description="Maximum number of results per page")
    page_no: Optional[int] = Field(1, description="Page number to retrieve (starting from 1)")
    filter_min: Optional[float] = Field(0, description="Minimum max_het_freq (inclusive) filter applied after caching")
    filter_max: Optional[float] = Field(1, description="Maximum max_het_freq (inclusive) filter applied after caching")

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
    page_limit: Optional[int] = Field(20, description="Maximum number of results per page")
    page_no: Optional[int] = Field(1, description="Page number to retrieve (starting from 1)")
    filter_min: Optional[float] = Field(0, description="Minimum pair_het_freq (inclusive) filter applied after caching")
    filter_max: Optional[float] = Field(1, description="Maximum pair_het_freq (inclusive) filter applied after caching")

class RegionResponse(BaseModel):
    """Response model for candidate regions"""
    results: List[RegionRecord]
    total_count: int
    page_no: int
    page_limit: int
    total_pages: int
    missing_transcripts: Optional[List[str]] = None
    exon_table: List[ExonRecord]

class HetFreqResponse(BaseModel):
    """Response model for heterozygous frequency results"""
    results: List[HetFreqRecord]
    total_count: int
    page_no: int
    page_limit: int
    total_pages: int
    population: str

class PairHetFreqResponse(BaseModel):
    """Response model for pair heterozygous frequency results"""
    results: List[PairHetFreqRecord]
    total_count: int
    page_no: int
    page_limit: int
    total_pages: int

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

# Helper function to get exons for transcript IDs
async def get_exons_for_transcripts(transcript_ids: List[str], conn: Optional[sqlite3.Connection] = None) -> List[ExonRecord]:
    """Get exon information for the given transcript IDs"""
    should_close_conn = conn is None
    if conn is None:
        conn = get_db_connection()

    try:
        # Query all exons for the given transcript_ids with transcript_id mapping
        if transcript_ids:
            # Build placeholders for the IN clause
            placeholders = ','.join(['?' for _ in transcript_ids])

            # First get transcript to exon mapping
            transcript_query = f"""
            SELECT transcript_id, exon_list, strand
            FROM transcripts
            WHERE transcript_id IN ({placeholders})
            """
            df_transcripts = pd.read_sql(transcript_query, conn, params=transcript_ids)

            exon_ids = []
            transcript_to_exons = pd.DataFrame()

            for _, row in df_transcripts.iterrows():
                transcript_id = row['transcript_id']
                if pd.notna(row['exon_list']) and row['exon_list']:
                    exons = row['exon_list'].split(',')
                    exon_ids.extend(exons)
                    transcript_to_exons = pd.concat([transcript_to_exons,
                                                     pd.DataFrame({'transcript_id': transcript_id,
                                                                   'exon_id': exons})],
                                                    ignore_index=True)

            if exon_ids:
                exon_placeholders = ','.join(['?' for _ in exon_ids])
                exon_query = f"""
                SELECT exon_id, start, end
                FROM exons
                WHERE exon_id IN ({exon_placeholders})
                """
                df_exons = pd.read_sql(exon_query, conn, params=exon_ids)

                # Use pandas join to efficiently map transcript_id to exons
                df_exons = df_exons.merge(transcript_to_exons, on='exon_id', how='left')

                # annotate exon number within each
                # if strand is +, count by increasing start position
                # if strand is -, count by decreasing start position
                tx_strand = df_transcripts['strand'].iloc[0]
                if tx_strand == '+':
                    df_exons['exon_no'] = df_exons.groupby('transcript_id')['start'].rank(method='dense', ascending=True)
                elif tx_strand == '-':
                    df_exons['exon_no'] = df_exons.groupby('transcript_id')['start'].rank(method='dense', ascending=False)
                else:
                    logger.error(f"Unknown transcript strand: {tx_strand}")

            else:
                df_exons = pd.DataFrame()
        else:
            df_exons = pd.DataFrame()

        # sort exon by transcript_id and exon_no
        df_exons = df_exons.sort_values(by=['transcript_id', 'exon_no'])

        # Convert to response models
        exon_table = [ExonRecord(**row) for row in df_exons.to_dict('records')]
        return exon_table

    finally:
        if should_close_conn and conn:
            conn.close()
            
@app.get("/")
async def root():
    """Root endpoint with API information"""
    logger.info("Root endpoint called - returning API information")
    return {
        "message": "COVER API",
        "version": "0.1.0",
        "endpoints": {
            "/query_results": "Query het_freq results by gene_id, transcript_id, or variant_id (ResultsRequest/ResultsResponse)",
            "/step1": "Query database by gene_id or gene_name (TranscriptRequest/TranscriptResponse)",
            "/step2": "Find candidate regions for transcript IDs with exon table (RegionRequest/RegionResponse)",
            "/step3": "Calculate heterozygous frequencies for SNP pairs (HetFreqRequest/HetFreqResponse)",
            "/step4": "Calculate pair heterozygous frequencies for combinations (PairHetFreqRequest/PairHetFreqResponse)",
            "/health": "Health check endpoint"
        }
    }

@app.get("/health")
async def health_check():
    """Health check endpoint"""
    logger.info("Health check endpoint called - checking system status")
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
        ResultsResponse containing filtered heterozygous frequency results with pagination
    """
    logger.info(f"Query results - gene_id: {request.gene_id}, transcript_id: {request.transcript_id}, variant_id: {request.variant_id}")
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

        # Initialize pagination parameters
        page_limit, page_no = initialize_pagination_params(request.page_limit, request.page_no)

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

        # Get total count for pagination
        count_query = f"""
        SELECT COUNT(*) as total_count
        FROM het_freq
        WHERE {where_clause}
        """
        df_count = pd.read_sql(count_query, conn, params=params)
        total_count = df_count.iloc[0]['total_count'] if not df_count.empty else 0

        # Calculate pagination values
        total_pages = (total_count + page_limit - 1) // page_limit  # Ceiling division
        offset = (page_no - 1) * page_limit

        # Query het_freq data with pagination
        query = f"""
        SELECT transcript_id, gene_id, gene_name, variant1, variant1_region,
               variant2, variant2_region, distance, target, consequence,
               cis_het_freq, trans_het_freq, max_het_freq, target_genotype, population
        FROM het_freq
        WHERE {where_clause}
        LIMIT ? OFFSET ?
        """
        params.append(page_limit)
        params.append(offset)

        df_results = pd.read_sql(query, conn, params=params)

        # Convert to response models
        if not df_results.empty:
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
                'target_genotype': 'string',
                'population': 'string'
            })

            # sort by max_het_freq in descending order
            df_results = df_results.sort_values(by='max_het_freq', ascending=False)

            # Convert to list of dictionaries using pandas' optimized method
            results = [HetFreqRecord(**row) for row in df_results.to_dict('records')]
        else:
            results = []

        return ResultsResponse(
            results=results,
            total_count=total_count,
            page_no=page_no,
            page_limit=page_limit,
            total_pages=total_pages,
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
    logger.info(f"Step1 - gene_id: {request.gene_id}, gene_name: {request.gene_name}")
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
        SELECT ROW_NUMBER() OVER () as transcript_no,
               t.transcript_id, t.gene_id, t.gene_name, t.exon_list,
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
            # Extract transcript_ids for exon querying
            transcript_ids = df_transcripts['transcript_id'].tolist()

            # Use the shared function to get exons
            exon_table = await get_exons_for_transcripts(transcript_ids, conn)
        else:
            exon_table = []

        # Remove exon_list from transcript data before returning
        if not df_transcripts.empty:
            df_transcripts = df_transcripts.drop(columns=['exon_list'])

        # Convert to response models
        transcript_table = [TranscriptRecord(**row) for row in df_transcripts.to_dict('records')]
        # exon_table is already a list of ExonRecord objects from get_exons_for_transcripts()

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
    Find candidate regions for given transcript_id(s) with pagination support and exon table.

    Results are computed once and cached for subsequent pagination requests.
    Also returns exon information for the corresponding transcripts.

    Example usage:
        "transcript_ids": ["ENST00000354192", "ENST00000348990"]

    Returns:
        RegionResponse containing paginated candidate regions and exon table as JSON
    """
    logger.info(f"Step2 - transcript_ids: {request.transcript_ids}")
    try:
        # Normalize transcript_ids to list
        if isinstance(request.transcript_ids, str):
            transcript_ids = [request.transcript_ids]
        else:
            transcript_ids = request.transcript_ids

        if not transcript_ids:
            raise HTTPException(status_code=400, detail="At least one transcript_id must be provided")

        # Initialize pagination parameters
        page_limit, page_no = initialize_pagination_params(request.page_limit, request.page_no)

        # Increment cleanup counter for cache management
        increment_cleanup_counter()

        # Create cache key based on transcript_ids and parameters
        cache_params = {
            'transcript_ids': sorted(transcript_ids),
            'max_deletion': request.max_deletion,
            'splice_donor_len': request.splice_donor_len,
            'splice_receptor_len': request.splice_receptor_len,
            'n_before_stop': request.n_before_stop
        }
        cache_key = create_cache_key(cache_params, 'step2')
        cache_path = get_cache_path(cache_key, 'step2')

        logger.debug(f"Cache key: {cache_key}")
        logger.debug(f"Cache path: {cache_path}")
        logger.debug(f"Cache valid: {is_cache_valid(cache_path)}")

        # Save input parameters to cache directory
        save_input_parameters(request, cache_path)

        # Check if results are already cached
        if is_cache_valid(cache_path, ['candidate_regions.candidate_region.txt']):
            logger.info(f"Using cached results for cache key: {cache_key}")

            cached_results = load_cached_results(cache_path, {
                'results': 'candidate_regions.candidate_region.txt',
                'missing_transcripts': 'candidate_regions.missing_transcript.txt'
            })

            all_results = cached_results.get('results', [])
            missing_transcripts = cached_results.get('missing_transcripts')

            # Query exons for the transcript_ids
            exon_table = await get_exons_for_transcripts(transcript_ids, conn=None)

            # Use the common pagination function (handles defaults and validation)
            pagination_info = get_pagination_info(all_results, page_limit, page_no)

            logger.debug(f"Returning {len(pagination_info['results'])} results for page {page_no}")

            return RegionResponse(
                results=pagination_info['results'],
                total_count=pagination_info['total_count'],
                page_no=pagination_info['page_no'],
                page_limit=pagination_info['page_limit'],
                total_pages=pagination_info['total_pages'],
                missing_transcripts=missing_transcripts,
                exon_table=exon_table
            )

        # Results not cached, need to compute them
        logger.info(f"Computing candidate regions for cache key: {cache_key}")

        # Create temporary file for input
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as temp_input:
            temp_input.write('\n'.join(transcript_ids))
            temp_input_path = temp_input.name

        try:
            # Write directly to cache directory
            output_prefix = os.path.join(cache_path, 'candidate_regions')

            # Create arguments namespace for main_find_candidate_region
            args = argparse.Namespace(
                id=None,  # Using id_list instead
                id_list=temp_input_path,
                db=TX_DATABASE_PATH,
                output=output_prefix,
                max_deletion=request.max_deletion,
                splice_donor_len=request.splice_donor_len,
                splice_receptor_len=request.splice_receptor_len,
                n_before_stop=request.n_before_stop,
                include_start_loss=True, # for API server, alwasy include start loss
                start_loss_only=False
            )

            # Run the main function
            main_find_candidate_region(args)

            # Read the output files directly from cache
            candidate_file = f"{output_prefix}.candidate_region.txt"
            missing_file = f"{output_prefix}.missing_transcript.txt"

            all_results = []
            missing_transcripts = None

            # Read candidate regions if file exists
            if os.path.exists(candidate_file):
                df_results = pd.read_csv(candidate_file, sep='\t')
                all_results = [RegionRecord(**row) for row in df_results.to_dict('records')]
            else:
                logger.warning(f"Candidate file not found: {candidate_file}")

            # Read missing transcripts if file exists
            if os.path.exists(missing_file):
                with open(missing_file, 'r') as f:
                    content = f.read().strip()
                    if content:
                        missing_transcripts = [line.strip() for line in content.split('\n') if line.strip()]
                    else:
                        missing_transcripts = []  # Empty list for empty file
            else:
                logger.debug(f"Missing file not found: {missing_file} (this is normal)")

        finally:
            # Clean up temporary input file
            os.unlink(temp_input_path)

        # Query exons for the transcript_ids
        exon_table = await get_exons_for_transcripts(transcript_ids, conn=None)

        # Use the common pagination function
        pagination_info = get_pagination_info(all_results, page_limit, page_no)

        logger.debug(f"Returning {len(pagination_info['results'])} results for page {page_no} from computation")

        return RegionResponse(
            results=pagination_info['results'],
            total_count=pagination_info['total_count'],
            page_no=pagination_info['page_no'],
            page_limit=pagination_info['page_limit'],
            total_pages=pagination_info['total_pages'],
            missing_transcripts=missing_transcripts,
            exon_table=exon_table
        )

    except Exception as e:
        logger.error(f"Error in get_candidate_regions: {str(e)}")
        raise HTTPException(status_code=500, detail=f"Region analysis failed: {str(e)}")

@app.post("/step3", response_model=HetFreqResponse)
async def get_heterozygous_frequencies(request: HetFreqRequest):
    """
    Calculate heterozygous frequencies for SNP pairs in candidate regions with pagination support.

    Results are computed once and cached for subsequent pagination requests.
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
        HetFreqResponse containing paginated heterozygous frequency results as JSON
    """
    logger.info(f"Step3 - number of regions: {len(request.regions)}, population: {request.population}, " + 
                f"max_deletion: {request.max_deletion}, maf: {request.maf}, exc_het: {request.exc_het}, n_pair_max: {request.n_pair_max}")
    # Validate VCF file is configured
    if VCF_FILE_PATH is None:
        raise HTTPException(status_code=500, detail="VCF file not configured")
    if not os.path.exists(VCF_FILE_PATH):
        raise HTTPException(status_code=500, detail=f"VCF file not found: {VCF_FILE_PATH}")

    try:
        # Validate input
        if not request.regions:
            raise HTTPException(status_code=400, detail="At least one region must be provided")

        # Initialize pagination parameters
        page_limit, page_no = initialize_pagination_params(request.page_limit, request.page_no)

        # Increment cleanup counter for cache management
        increment_cleanup_counter()

        # Create cache key based on regions and parameters
        cache_params = {
            'regions': sorted([str(region) for region in request.regions]),  # Convert regions to strings for hashing
            'population': request.population,
            'max_deletion': request.max_deletion,
            'maf': request.maf,
            'exc_het': request.exc_het,
            'n_pair_max': request.n_pair_max
        }
        cache_key = create_cache_key(cache_params, 'step3')
        cache_path = get_cache_path(cache_key, 'step3')

        logger.debug(f"Cache key: {cache_key}")
        logger.debug(f"Cache path: {cache_path}")
        logger.debug(f"Cache valid: {is_cache_valid(cache_path)}")

        # Save input parameters to cache directory
        save_input_parameters(request, cache_path)

        # Check if results are already cached
        if is_cache_valid(cache_path, ['het_freq.het_freq.txt']):
            logger.info(f"Using cached results for cache key: {cache_key}")

            cached_results = load_cached_results(cache_path, {
                'results': 'het_freq.het_freq.txt'
            })

            # Get cached results (DataFrame expected)
            df_results = cached_results.get('results', None)
            df_results = add_het_freq_meta(df_results, request.regions, request.population)
            logger.debug(f"Cached results count: {len(df_results) if isinstance(df_results, pd.DataFrame) else 0}")
        # Results not cached, run het freq calculation
        else:
            df_results = run_het_freq_calculation(
                regions=request.regions,
                population=request.population,
                cache_key=cache_key,
                cache_path=cache_path,
                pair_het_cutoff=1.0,  # Force to 1 to skip pair calculations
                top_n_comb=1000,
                output_file_suffix="het_freq.txt",
                process_pair_results=False
            )
            
            if df_results.empty:
                raise HTTPException(status_code=500, detail="No available SNPs targeting the selected regions, please expand your selection.")
                        
        # Apply filters on max_het_freq
        if not df_results.empty:
            mask = (df_results['max_het_freq'].astype(float) >= request.filter_min) & \
                    (df_results['max_het_freq'].astype(float) <= request.filter_max)
            df_results = df_results[mask]
        
        # sort by max_het_freq in descending order
        df_results = df_results.sort_values(by='max_het_freq', ascending=False)

        # Use the common pagination function with DataFrame
        pagination_info = get_pagination_info(df_results, page_limit, page_no)

        logger.debug(f"Returning {len(pagination_info['results'])} results for page {page_no}")
        return HetFreqResponse(
            results=pagination_info['results'],
            total_count=pagination_info['total_count'],
            page_no=pagination_info['page_no'],
            page_limit=pagination_info['page_limit'],
            total_pages=pagination_info['total_pages'],
            population=request.population
        )

    except Exception as e:
        # Re-raise HTTPExceptions
        if isinstance(e, HTTPException):
            raise e
        
        logger.error(f"Error in get_heterozygous_frequencies: {str(e)}")
        raise HTTPException(status_code=500, detail=f"Heterozygous frequency calculation failed: {str(e)}")

@app.post("/step4", response_model=PairHetFreqResponse)
async def get_pair_heterozygous_frequencies(request: PairHetFreqRequest):
    """
    Calculate pair heterozygous frequencies for combinations of SNP pairs in candidate regions.
    This function allows user to specify pair_het_cutoff to control which pairs are considered.
    Results can be filtered by pair_het_freq using filter_min and filter_max parameters.

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
    logger.info(f"Step4 - number of regions: {len(request.regions)}, pair_het_cutoff: {request.pair_het_cutoff}" + 
                f", n_pair_max: {request.n_pair_max}, top_n_comb: {request.top_n_comb}" +
                f", filter_min: {request.filter_min}, filter_max: {request.filter_max}")
    # Validate VCF file is configured
    if VCF_FILE_PATH is None:
        raise HTTPException(status_code=500, detail="VCF file not configured")
    if not os.path.exists(VCF_FILE_PATH):
        raise HTTPException(status_code=500, detail=f"VCF file not found: {VCF_FILE_PATH}")

    try:
        # Validate input
        if not request.regions:
            raise HTTPException(status_code=400, detail="At least one region must be provided")

        # Initialize pagination parameters
        page_limit, page_no = initialize_pagination_params(request.page_limit, request.page_no)

        # Increment cleanup counter for cache management
        increment_cleanup_counter()

        # Create cache key based on regions and parameters
        cache_params = {
            'regions': sorted([str(region) for region in request.regions]),  # Convert regions to strings for hashing
            'population': request.population,
            'pair_het_cutoff': request.pair_het_cutoff,
            'max_deletion': request.max_deletion,
            'maf': request.maf,
            'exc_het': request.exc_het,
            'n_pair_max': request.n_pair_max,
            'top_n_comb': request.top_n_comb
        }
        cache_key = create_cache_key(cache_params, 'step4')
        cache_path = get_cache_path(cache_key, 'step4')

        logger.debug(f"Cache key: {cache_key}")
        logger.debug(f"Cache path: {cache_path}")
        logger.debug(f"Cache valid: {is_cache_valid(cache_path)}")

        # Save input parameters to cache directory
        save_input_parameters(request, cache_path)

        # Check if results are already cached
        if is_cache_valid(cache_path, ['pair_het_freq.pair_het_freq.txt']):
            logger.info(f"Using cached results for cache key: {cache_key}")

            cached_results = load_cached_results(cache_path, {
                'results': 'pair_het_freq.pair_het_freq.txt'
            })

            # Get cached results
            df_results = cached_results.get('results', None)
            df_results = add_het_freq_meta(df_results, request.regions, request.population)
            logger.debug(f"Cached results count: {len(df_results)}")
        # Results not cached, run het freq calculation
        else:
            df_results = run_het_freq_calculation(
                regions=request.regions,
                population=request.population,
                cache_key=cache_key,
                cache_path=cache_path,
                pair_het_cutoff=request.pair_het_cutoff,
                top_n_comb=request.top_n_comb,
                output_file_suffix="pair_het_freq.txt",
                process_pair_results=True
            )

        # Apply filters on pair_het_freq
        if not df_results.empty:
            mask = (df_results['pair_het_freq'].astype(float) >= request.filter_min) & \
                    (df_results['pair_het_freq'].astype(float) <= request.filter_max)
            df_results = df_results[mask]
        
        # sort by pair_het_freq in descending order
        df_results = df_results.sort_values(by='pair_het_freq', ascending=False)

        # Use the common pagination function
        pagination_info = get_pagination_info(df_results, page_limit, page_no)

        logger.debug(f"Returning {len(pagination_info['results'])} results for page {page_no}")

        return PairHetFreqResponse(
            results=pagination_info['results'],
            total_count=pagination_info['total_count'],
            page_no=pagination_info['page_no'],
            page_limit=pagination_info['page_limit'],
            total_pages=pagination_info['total_pages']
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
    parser.add_argument("--allowed-origins", default=None, help="Comma separated list of allowed origins")
    parser.add_argument("--debug", action="store_true", help="Enable debug logging")

    args = parser.parse_args()

    # Configure logging level based on debug flag
    if args.debug:
        # Update existing logger configuration for debug mode
        logger.setLevel(logging.DEBUG)
        for handler in logger.handlers:
            handler.setLevel(logging.DEBUG)
        logger.info("Debug logging enabled")

    # Log the current log file location
    current_date = datetime.datetime.now().strftime('%Y-%m-%d')
    log_filename = f'cover-api.{current_date}.log'
    logs_dir = os.path.join(os.getcwd(), 'logs')
    log_filepath = os.path.join(logs_dir, log_filename)
    logger.info(f"Log file: {log_filepath}")

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

    # Configure CORS
    if args.allowed_origins:
        lst_allowed_origins = args.allowed_origins.split(",")
        configure_cors(app, lst_allowed_origins)
        logger.info(f"Allowed origins: {lst_allowed_origins}")

    uvicorn.run(app, host=args.host, port=args.port)

if __name__ == "__main__":
    main()