#!/usr/bin/env python
"""
Convert a GenBank reference into the FASTA and GFF3 files for downstream workflow steps e,g. makeblastdb, indexing.

Last checked by: BP
"""
from Bio import SeqIO
from Bio.Seq import UndefinedSequenceError
import argparse
import sys
import multiprocessing as mp

def get_strand_symbol(strand):
    """Get strand symbols."""
    return '+' if strand == 1 else '-' if strand == -1 else '.'

def get_frame(feature):
    """Get frame for CDS features."""
    if feature.type == 'CDS':
        return str(feature.location.start % 3)
    return '.'

def format_attributes(feature):
    """Build the GFF attributes."""
    attributes = []
    
    #ID is required - construct from locus_tag or gene
    feature_id = (
        feature.qualifiers.get('locus_tag', [''])[0] or
        feature.qualifiers.get('gene', [''])[0] or
        f"{feature.type}_{int(feature.location.start)+1}_{int(feature.location.end)}"
    )
    attributes.append(f"ID={feature_id}")
    
    #add other common qualifiers
    if 'gene' in feature.qualifiers:
        attributes.append(f"Gene={feature.qualifiers['gene'][0]}")
    if 'product' in feature.qualifiers:
        attributes.append(f"Product={feature.qualifiers['product'][0]}")
    
    return ';'.join(attributes)

def process_record(record):
    """Extract sequence and attributes from selected GenBank feature record for ref. track annotation."""
    included_features = ["CDS", "tRNA", "rRNA", "ncRNA", "misc_feature", "mobile_element", "source"]
    record_features = []
    
    try:
        # Verify sequence is defined
        str(record.seq)
        print(f"Processing record {record.id}")
        
        # Extract features
        for feature in record.features:
            if feature.type in included_features:
                try:
                    feature_seq = feature.extract(record.seq)
                    record_features.append({
                        'id': f"{record.id}_{int(feature.location.start)+1}_{int(feature.location.end)}",
                        'seq': feature_seq,
                        'feature': feature,
                        'record_id': record.id
                    })
                except Exception as e:
                    print(f"Warning: Error processing feature in record {record.id}: {e}")
                    continue
        
        return {'id': record.id, 'seq': record.seq, 'features': record_features}
        
    except UndefinedSequenceError as e:
        print(f"Error (UndefinedSequenceError): {e} in record {record.id}")
        return None
    except Exception as e:
        print(f"Error: {e} in record {record.id}")
        return None

def setup_reference(ref_gbk, out_seqs, out_features, threads='max'):
    """Convert GenBank reference input into FASTA and features outputs."""
    # Determine thread count
    if threads == 'max' or threads is None:
        n_cpus = max(1, mp.cpu_count() - 2)  # Reserve 2 CPUs for overhead
    else:
        try:
            n_cpus = int(threads)
            n_cpus = max(1, min(n_cpus, mp.cpu_count()))
        except ValueError:
            print(f"Warning: Invalid threads value '{threads}', using maximum available minus 2")
            n_cpus = max(1, mp.cpu_count() - 2)
            
    print(f"Using {n_cpus} CPU cores")
    
    try:
        # Parse all records (must be done sequentially)
        records = list(SeqIO.parse(ref_gbk, "genbank"))
        
        if not records:
            print("Error: No records found in the GenBank file")
            sys.exit(1)
            
        print(f"Found {len(records)} records in GenBank file")
        
        # Process records in parallel
        ref_feature_dict = {}
        ref_seqs = {}
        
        with mp.Pool(processes=n_cpus) as pool:
            results = pool.map(process_record, records)
            
            # Collect results
            for result in results:
                if result is None:
                    continue
                    
                if result['id'] in ref_seqs:
                    print(f"Warning: Sequence {result['id']} is already present in the reference genome, will use the first one")
                else:
                    ref_seqs[result['id']] = result['seq']
                    
                # Add features
                for feature_data in result['features']:
                    ref_feature_dict[feature_data['id']] = {
                        'seq': feature_data['seq'],
                        'feature': feature_data['feature'],
                        'record_id': feature_data['record_id']
                    }
        
        if not ref_feature_dict:
            print("Error: No features found in the reference genome")
            sys.exit(1)

        print(f"Writing {len(ref_seqs)} sequences to {out_seqs}")
        # Write out sequences
        with open(out_seqs, "w") as f:
            for rec_id, seq in ref_seqs.items():
                f.write(f">{rec_id}\n{seq}\n")

        print(f"Writing {len(ref_feature_dict)} features to {out_features}")
        # Write out GFF3 file
        with open(out_features, "w") as f:
            # Header
            f.write("##gff-version 3\n")
            
            # Features
            for feature_id, data in ref_feature_dict.items():
                feature = data['feature']
                record_id = data['record_id']
                fields = [
                    record_id,                                   # 1. seqid
                    "GenBank",                                   # 2. source
                    feature.type,                                # 3. type
                    str(int(feature.location.start) + 1),        # 4. start (1-based)
                    str(int(feature.location.end)),              # 5. end (1-based)
                    ".",                                         # 6. score
                    get_strand_symbol(feature.location.strand),  # 7. strand
                    get_frame(feature),                          # 8. frame
                    format_attributes(feature)                   # 9. attributes
                ]
                
                f.write('\t'.join(fields) + '\n')

        print("Conversion completed successfully")
        
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    """Parse args and generate reference FASTA and GFF3 files."""
    parser = argparse.ArgumentParser(description="Convert GenBank to FASTA and GFF3 files with parallel processing")
    parser.add_argument("-r", "--ref_gbk", help="A GenBank reference file (.gb, .gbk)", required=True)
    parser.add_argument("-s", "--out_seqs", help="Output multi-FASTA file for sequences", default="ref.fna")
    parser.add_argument("-f", "--out_features", help="Output GFF3 file with feature annotations", default="ref_features.gff3")
    parser.add_argument("-t", "--threads", help="Number of CPU threads to use (default: max, using all CPUs minus 2)", default="max")
    args = parser.parse_args()
    
    setup_reference(args.ref_gbk, args.out_seqs, args.out_features, args.threads)
