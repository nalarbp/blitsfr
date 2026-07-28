#!/usr/bin/env python
"""
Build CGView-compatible JSON for BLITSFR reports.

This is Python-based reimplementation of the map-building workflow used by CGView.js. 
It prepares CGView-compatible JSON for the bundled viewer.

Upstream project references:
- CGView.js docs: https://js.cgview.ca/
- CGView.js repository: https://github.com/sciguy/cgview-js

Last checked by: BP
"""
import argparse
import json
import uuid
import sys
import re
import time
import psutil
import os
from collections import defaultdict
from datetime import datetime
from dataclasses import dataclass, asdict, field
from typing import List, Dict, Union, Optional, Tuple, Any
from pathlib import Path

#Todo: to_dict has the highest memory and time usage (tried, keep it for this version, rebase for next milestone)

def track_performance(func):
    """Used for perf. eval, getting timing and memory logs."""
    def wrapper(*args, **kwargs):
        # Large report builds can be memory-heavy, so keep lightweight timing telemetry here.
        process = psutil.Process(os.getpid())
        start_memory = process.memory_info().rss / 1024 / 1024  # MB
        start_time = time.time()
        
        result = func(*args, **kwargs)
        
        end_time = time.time()
        end_memory = process.memory_info().rss / 1024 / 1024  # MB
        elapsed_time = end_time - start_time
        memory_used = end_memory - start_memory
        
        print(f"Function {func.__name__} completed in {elapsed_time:.4f} seconds")
        print(f"Memory change: {memory_used:.2f} MB, Total memory: {end_memory:.2f} MB")
        
        return result
    return wrapper

@dataclass
class Settings:
    format: str = "circular"
    geneticCode: int = 11
    backgroundColor: str = "rgba(255,255,255,1)"
    showShading: bool = True
    arrowHeadLength: float = 0.3
    initialMapThicknessProportion: float = 0.25
    maxMapThicknessProportion: float = 0.4
            
    def to_dict(self):
        return asdict(self)

@dataclass
class Backbone:
    color: str = "rgba(128,128,128,1)"
    colorAlternate: str = "rgba(200,200,200,1)"
    thickness: int = 5
    decoration: str = "arrow"
    visible: bool = True

    def to_dict(self):
        return asdict(self)

@dataclass
class Ruler:
    font: str = "sans-serif,plain,10"
    color: str = "rgba(0,0,0,1)"
    visible: bool = True

    def to_dict(self):
        return asdict(self)

@dataclass
class Annotation:
    font: str = "monospace,plain,12"
    color: str = "black"
    onlyDrawFavorites: bool = False
    labelPlacement: str = "default"
    visible: bool = True

    def to_dict(self):
        return asdict(self)

@dataclass
class Divider:
    visible: bool = False
    color: str = "rgba(0,0,0,1)"
    thickness: int = 0
    spacing: int = 0

    def to_dict(self):
        return asdict(self)

@dataclass
class Dividers:
    slot: Divider = field(default_factory=Divider)
    track: Divider = field(default_factory=Divider)

    def to_dict(self):
        return {
            'slot': self.slot.to_dict(),
            'track': self.track.to_dict()
        }
        
@dataclass
class Highlighter:
    visible: bool = True

    def to_dict(self):
        return asdict(self)

@dataclass
class Captions:
    name: str = "CGView map"
    position: str = "bottom-center"
    textAlignment: str = "center"
    font: str = "sans-serif,plain,24"
    fontColor: str = "rgba(0,0,139,1)"
    backgroundColor: str = "rgba(255,255,255,1)"

    def to_dict(self):
        return asdict(self)

@dataclass
class LegendItem:
    name: str
    swatchColor: str
    decoration: str = "arc"
    
    def to_dict(self):
        return asdict(self)
    
@dataclass
class Legend:
    position: Union[str, dict] = "top-right"
    anchor: Union[str, dict] = "auto"
    textAlignment: str = "left"
    defaultFont: str = "sans-serif,plain,14"
    defaultFontColor: str = "rgba(0,0,0,1)"
    defaultMinArcLength: float = 1
    backgroundColor: str = "rgba(255,255,255,0.75)"
    visible: bool = True
    meta: dict = field(default_factory=dict)
    items: List[LegendItem] = field(default_factory=list)
    
    def to_dict(self):
        return {
            **{k: v for k, v in asdict(self).items() if k != 'items'},
            'items': [item.to_dict() for item in self.items]
        }

@dataclass
class Contig:
    name: str
    orientation: str  # forward or reverse
    length: int
    seq: str = ""
    visible: bool = True

    def to_dict(self):
        return asdict(self)

@dataclass
class Sequence:
    font: str = "sans-serif,plain,14"
    color: str = "rgb(0,0,0)"
    contigs: List[Contig] = field(default_factory=list)
    
    def add_contig(self, contig: Contig):
        """Append a contig to sequence."""
        self.contigs.append(contig)
    
    def to_dict(self):
        return {
            **{k: v for k, v in asdict(self).items() if k != 'contigs'},
            'contigs': [contig.to_dict() for contig in self.contigs]
        }
    
@dataclass
class Feature:
    name: str
    type: str
    legend: Union[str, dict]
    source: str
    contig: str 
    start: int
    stop: int
    strand: str
    tags: List[str] = field(default_factory=list)
    score: Optional[float] = None
    favorite: bool = False
    visible: bool = True
    meta: dict = field(default_factory=dict)

    def __post_init__(self):
        """Check and normalize feature coordinate order."""
        # Normalize reversed coordinates so downstream rendering always sees start <= stop.
        if self.start > self.stop:
            self.start, self.stop = self.stop, self.start
        if self.strand not in ['+', '-', '.', 1, -1]:
            raise ValueError(f"Invalid strand value: {self.strand}")

    def to_dict(self):
        return {k: v for k, v in asdict(self).items() if v is not None}

@dataclass
class Plot:
    name: str
    positions: List[int]
    scores: List[float]
    source: str
    baseline: float = 0
    legend: str = ""
    axisMin: float = 0
    visible: bool = True
    
    def __post_init__(self):
        """Validate plot array lengths and default the legend label."""
        if len(self.positions) != len(self.scores):
            raise ValueError("Number of positions must match number of scores")
        if not self.legend:
            self.legend = self.name
    
    def to_dict(self):
        return asdict(self)
    
@dataclass
class Track:
    name: str
    dataType: str
    dataMethod: str
    dataKeys: Union[str, List[str]]
    position: str = "inside"
    separateFeaturesBy: str = "none"
    thicknessRatio: float = 1.0
    favorite: bool = False
    visible: bool = True
    meta: dict = field(default_factory=dict)

    def to_dict(self):
        return {k: v for k, v in asdict(self).items() if v is not None}

@dataclass
class CGViewBuilder:
    """Main CGView map instace."""
    name: str = 'CGViewMap'
    version: str = "1.5.1"
    settings: Settings = field(default_factory=Settings)
    backbone: Backbone = field(default_factory=Backbone)
    features: List[Feature] = field(default_factory=list)
    plots: List[Plot] = field(default_factory=list)
    tracks: List[Track] = field(default_factory=list)
    sequence: Optional[Sequence] = None
    legend: Legend = field(default_factory=Legend)
    annotation: Annotation = field(default_factory=Annotation)
    ruler: Ruler = field(default_factory=Ruler)
    captions: Captions = field(default_factory=Captions)
    dividers: Dividers = field(default_factory=Dividers)
    highlighter: Highlighter = field(default_factory=Highlighter)
    id: str = field(default_factory=lambda: uuid.uuid4().hex[:20])
    created: str = field(default_factory=lambda: datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    _feature_sources: set = field(default_factory=set)
    
    def __post_init__(self):
        self.captions.name = self.name

    def set_sequence(self, sequence: Sequence):
        """Add ref. sequence to the builder."""
        self.sequence = sequence

    def add_feature(self, feature: Feature):
        """Add a feature to the buildert."""
        self.features.append(feature)
        if feature.source == "genbank-features" and feature.source not in self._feature_sources:
            self._feature_sources.add(feature.source)
            track = Track(
                name="Features",
                dataType="feature",
                dataMethod="source",
                dataKeys=feature.source,
                position="both",
                separateFeaturesBy="strand",
                thicknessRatio=1
            )
            self.tracks.append(track)

    def add_plot(self, plot: Plot):
        """Add plot to the builder."""
        self.plots.append(plot)
        track = Track(
            name=plot.name,
            dataType="plot",
            dataMethod="source",
            dataKeys=plot.source,
            separateFeaturesBy="none",
            thicknessRatio=1
        )
        self.tracks.append(track)

    def add_track(self, track: Track):
        """Append track to the builder."""
        self.tracks.append(track)

    @track_performance
    def to_dict(self) -> dict:
        """Make the final CGView.js compatible JSON payload."""
        result = {
            "cgview": {
                "name": self.name,
                "version": self.version,
                "created": self.created,
                "id": self.id,
                "settings": self.settings.to_dict(),
                "backbone": self.backbone.to_dict(),
                "ruler": self.ruler.to_dict(),
                "annotation": self.annotation.to_dict(),
                "dividers": self.dividers.to_dict(),
                "highlighter": self.highlighter.to_dict(),
                "captions": self.captions.to_dict(),
                "legend": self.legend.to_dict(),
                "features": [f.to_dict() for f in self.features],
                "tracks": [t.to_dict() for t in self.tracks],
                "plots": [p.to_dict() for p in self.plots]
            }
        }
        
        if self.sequence:
            result["cgview"]["sequence"] = self.sequence.to_dict()
            
        return result

@track_performance
def parse_sequence_file(file_path: str) -> Sequence:
    """Parse a FASTA file into a CGView sequence, works on one or more contigs."""
    try:
        with open(file_path, 'r') as f:
            lines = f.readlines()
            if not lines:
                raise ValueError("Empty sequence file")
            
            contigs = []
            current_contig = []
            current_name = ""
            
            for line in lines:
                line = line.strip()
                if line.startswith('>'):
                    if current_name and current_contig:
                        seq = ''.join(current_contig)
                        contigs.append(Contig(
                            name=current_name,
                            orientation="+",
                            length=len(seq),
                            seq=seq
                        ))
                    current_name = line[1:].split()[0]
                    current_contig = []
                else:
                    current_contig.append(line)
            
            if current_name and current_contig:
                seq = ''.join(current_contig)
                contigs.append(Contig(
                    name=current_name,
                    orientation="+",
                    length=len(seq),
                    seq=seq
                ))
            
            return Sequence(contigs=contigs)
            
    except FileNotFoundError:
        print(f"Error: File '{file_path}' not found.")
        sys.exit(1)
        
def parse_gff_attributes(attr_string: str, feature_id_attr: str = 'first') -> str:
    """Extract feature attr. from a GFF ref"""
    attributes = {}
    for attr in attr_string.split(';'):
        if not attr.strip():
            continue
        if '=' in attr:
            key, value = attr.split('=', 1)
            attributes[key.strip()] = value.strip()
    
    if feature_id_attr == 'first':
        return next(iter(attributes.values())) if attributes else ''
    return attributes.get(feature_id_attr, next(iter(attributes.values())) if attributes else '')

def get_feature_color(feature_type: str) -> str:
    """Get the default color for ref. feature type."""
    color_map = {
        'CDS': 'rgba(255,0,0,1)',
        'gene': 'rgba(0,255,0,1)',
        'tRNA': 'rgba(0,0,255,1)',
        'rRNA': 'rgba(255,165,0,1)',
        'ncRNA': 'rgba(128,0,128,1)',
        'mobile_element': 'rgba(165,42,42,1)',
        'misc_feature': 'rgba(128,128,128,1)'
    }
    return color_map.get(feature_type, 'rgba(100,100,100,1)')

def get_feature_decoration(feature_type: str) -> str:
    """Get the legend decoration for feature type."""
    arrow_features = {'CDS', 'gene', 'tRNA', 'rRNA', 'ncRNA'}
    return 'arrow' if feature_type in arrow_features else 'arc'

@track_performance
def parse_feature_file(file_path: str, feature_types: set = None, feature_id_attr: str = 'first') -> tuple[List[Feature], List[LegendItem]]:
    """Parse GFF features and build legend items."""
    features = []
    seen_feature_types = set()
    legend_items = []
    
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
                
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
                
            feature_type = fields[2]
            if feature_types and feature_type not in feature_types:
                continue
                
            try:
                feature_id = parse_gff_attributes(fields[8], feature_id_attr)
                if not feature_id:
                    print(f"Warning: No valid ID found in line: {line.strip()}")
                    continue
                
                color = get_feature_color(feature_type)
                
                feature = Feature(
                    name=feature_id,
                    type=feature_type,
                    legend=feature_type,
                    source="genbank-features",
                    contig=fields[0],
                    start=int(fields[3]),
                    stop=int(fields[4]),
                    strand=fields[6],
                    meta={
                        "name": feature_id,
                        "gene": parse_gff_attributes(fields[8], 'Gene'),
                        "product": parse_gff_attributes(fields[8], 'Product'),
                        "color": color}
                )
                features.append(feature)
                
                if feature_type not in seen_feature_types:
                    seen_feature_types.add(feature_type)
                    legend_items.append(LegendItem(
                        name=feature_type,
                        swatchColor=color,
                        decoration=get_feature_decoration(feature_type)
                    ))
                    
            except (ValueError, IndexError) as e:
                print(f"Warning: Skipping malformed feature line: {line.strip()}")
                continue
    
    return features, legend_items

@track_performance
def parse_blast_file(file_path: str) -> Tuple[Dict[str, List[Feature]], Dict[str, str]]:
    """Convert filtered BLAST rows into CGView feature tracks grouped by query."""
    pattern = re.compile(r'[^\w-]')
    features_by_query = defaultdict(list)
    query_files = set()
    next_line_is_header = True
    
    with open(file_path, 'r') as f:
        for line in f:
            if next_line_is_header:
                next_line_is_header = False
                continue
                
            fields = line.strip().split('\t')
            if len(fields) < 17:
                continue
                
            query_files.add(fields[0])
            query_file = fields[0]
            sseqid = fields[2]
            pident = float(fields[3])
            start = int(fields[9])
            end = int(fields[10])
            
            if len(fields) <= 17:
                continue
                
            if start < end:
                strand = 1  # "+"
            else:
                strand = -1  # "-"
                start, end = end, start
                
            # Keep the canonical query ID so report payloads and CGView tracks use the same key.
            genome_name = query_file
            # Previous sanitised track naming kept here for reference in case CGView-specific issues return.
            # genome_name = pattern.sub('_', query_file)
            
            feature = Feature(
                name="",
                type="blast",
                legend=genome_name,
                source="",  # Temporary placeholder, will be updated later
                contig=sseqid,
                start=start,
                stop=end,
                strand=strand,
                meta={
                    "identity": pident,
                    "mimatches": 0,
                    "evalue": 0,
                    "bit_score": 0
                }
            )
            
            features_by_query[genome_name].append(feature)
    
    # Preserve the original query identity before converting it into CGView track/source keys.
    query_to_source = {
        qfile: f"blast-{i+1}.1"
        for i, qfile in enumerate(sorted(query_files))
    }
    
    for genome_name, features in features_by_query.items():
        query_file = next((qfile for qfile in query_files if qfile == genome_name), None)
        # Previous sanitised lookup kept here for reference.
        # query_file = next((qfile for qfile in query_files if pattern.sub('_', qfile) == genome_name), None)
        if query_file:
            source = query_to_source[query_file]
            for feature in features:
                feature.source = source
    
    return dict(features_by_query), query_to_source


@track_performance
def add_blast_tracks(builder: CGViewBuilder, blast_features: Dict[str, List[Feature]], query_to_source: Dict[str, str], colors: List[str] = None):
    """Attach BLAST-derived features, tracks, and legend items to the builder."""
    if colors is None:
        colors = ["rgba(0,0,0,1)"]
    
    arrow_features = {'CDS', 'gene', 'tRNA', 'rRNA', 'ncRNA'}
    for genome_name, features in blast_features.items():
        color = colors[0]
        
        for feature in features:
            builder.add_feature(feature)
        
        decoration = "arrow" if features[0].type in arrow_features else "arc"
        legend_item = LegendItem(
            name=genome_name,
            swatchColor=color,
            decoration=decoration
        )
        builder.legend.items.append(legend_item)
        
        source = features[0].source
        
        track = Track(
            name=genome_name,
            dataType="feature",
            dataMethod="source",
            dataKeys=source,
            position="inside",
            separateFeaturesBy="none",
            thicknessRatio=1
        )
        builder.tracks.append(track)

@track_performance
def parse_plot_file(file_path: str, sequence: Sequence = None) -> tuple[Dict[str, List[Plot]], Dict[str, str]]:
    """Parse windowed coverage data into CGView plot objects grouped by sample."""
    plots_by_sample = {}
    sample_files = set()
    
    try:
        with open(file_path, 'r') as f:
            header = f.readline().strip().split('\t')
            
            expected_columns = ['SampleID', 'Template', 'Window', 'Start', 'End', 'Value']
            if len(header) < 6 or not all(col in header for col in expected_columns):
                raise ValueError(f"Invalid header format. Expected columns: {expected_columns}")
            
            sample_idx = header.index('SampleID')
            template_idx = header.index('Template')
            window_idx = header.index('Window')
            start_idx = header.index('Start')
            end_idx = header.index('End')
            value_idx = header.index('Value')
            
            for line in f:
                if line.startswith('#'):
                    continue
                    
                fields = line.strip().split('\t')
                if len(fields) < 6:
                    continue
                    
                sample_id = fields[sample_idx]
                sample_files.add(sample_id)
            
            # Keep sample IDs stable while assigning compact internal source names for tracks.
            sample_to_source = {
                sample: f"plot-{i+1}" 
                for i, sample in enumerate(sorted(sample_files))
            }
            
            f.seek(0)
            next(f)
            
            contig_offsets = {}
            contig_lengths = {}
            contig_order = []
            if sequence:
                offset = 0
                for contig in sequence.contigs:
                    contig_offsets[contig.name] = offset
                    contig_lengths[contig.name] = contig.length
                    contig_order.append(contig.name)
                    offset += contig.length
            
            grouped_data = {}
            for line in f:
                if line.startswith('#'):
                    continue
                    
                fields = line.strip().split('\t')
                if len(fields) < 6:
                    continue
                    
                sample_id = fields[sample_idx]
                template = fields[template_idx]
                window = fields[window_idx]
                start = int(fields[start_idx])
                end = int(fields[end_idx])
                try:
                    value = float(fields[value_idx])
                except ValueError:
                    print(f"Warning: Skipping line with invalid value: {line.strip()}")
                    continue
                
                if sequence:
                    if template in contig_offsets:
                        offset = contig_offsets[template]
                        start += offset
                        end += offset
                    else:
                        # Fall back to local coordinates when plot templates do not match contig names.
                        print(f"Warning: Template '{template}' not found in sequence contigs. Using local coordinates.")
                
                if sample_id not in grouped_data:
                    grouped_data[sample_id] = {
                        'positions': [],
                        'scores': [],
                        'windows': [],
                        'ends': [],
                        'templates': []
                    }
                
                grouped_data[sample_id]['positions'].append(start)
                grouped_data[sample_id]['scores'].append(value)
                grouped_data[sample_id]['windows'].append(window)
                grouped_data[sample_id]['ends'].append(end)
                grouped_data[sample_id]['templates'].append(template)
            
            for sample_id, data in grouped_data.items():
                plot_name = sample_id
                source = sample_to_source[sample_id]
                
                if sample_id not in plots_by_sample:
                    plots_by_sample[sample_id] = []
                
                if sequence and len(contig_order) > 1:
                    # Rebuild points in contig order and insert boundary markers to avoid cross-contig lines.
                    template_data = {}
                    for i, template in enumerate(data['templates']):
                        if template not in template_data:
                            template_data[template] = {
                                'positions': [],
                                'scores': [],
                                'indices': []
                            }
                        template_data[template]['positions'].append(data['positions'][i])
                        template_data[template]['scores'].append(data['scores'][i])
                        template_data[template]['indices'].append(i)
                    
                    new_positions = []
                    new_scores = []
                    new_windows = []
                    new_ends = []
                    new_templates = []
                    
                    for template in contig_order:
                        if template not in template_data:
                            continue
                            
                        t_positions = template_data[template]['positions']
                        t_scores = template_data[template]['scores']
                        t_indices = template_data[template]['indices']
                        
                        for i, idx in enumerate(t_indices):
                            new_positions.append(data['positions'][idx])
                            new_scores.append(data['scores'][idx])
                            new_windows.append(data['windows'][idx])
                            new_ends.append(data['ends'][idx])
                            new_templates.append(template)
                        
                        if template != contig_order[-1]:
                            offset = contig_offsets[template]
                            length = contig_lengths[template]
                            end_pos = offset + length
                            
                            max_score = max(t_scores) if t_scores else 0
                            
                            new_positions.append(end_pos)
                            new_scores.append(max_score)
                            new_windows.append(f"boundary_{template}_end")
                            new_ends.append(end_pos)
                            new_templates.append(template)
                            
                            new_positions.append(end_pos)
                            new_scores.append(0)
                            new_windows.append(f"boundary_{template}_end_zero")
                            new_ends.append(end_pos)
                            new_templates.append(template)
                            
                            next_template = contig_order[contig_order.index(template) + 1]
                            next_offset = contig_offsets[next_template]
                            
                            new_positions.append(next_offset + 1)
                            new_scores.append(0)
                            new_windows.append(f"boundary_{next_template}_start_zero")
                            new_ends.append(next_offset + 1)
                            new_templates.append(next_template)
                            
                            if next_template in template_data and template_data[next_template]['positions']:
                                next_score = template_data[next_template]['scores'][0]
                            else:
                                next_score = 0
                                
                            new_positions.append(next_offset + 1)
                            new_scores.append(next_score)
                            new_windows.append(f"boundary_{next_template}_start")
                            new_ends.append(next_offset + 1)
                            new_templates.append(next_template)
                    
                    plot = Plot(
                        name=plot_name,
                        positions=new_positions,
                        scores=new_scores,
                        source=source,
                        legend=sample_id,
                        visible=True
                    )
                    
                    plot.meta = {
                        "windows": new_windows,
                        "ends": new_ends,
                        "templates": new_templates,
                        "has_boundaries": True
                    }
                else:
                    plot = Plot(
                        name=plot_name,
                        positions=data['positions'],
                        scores=data['scores'],
                        source=source,
                        legend=sample_id,
                        visible=True
                    )
                    
                    plot.meta = {
                        "windows": data['windows'],
                        "ends": data['ends'],
                        "templates": data['templates'],
                        "has_boundaries": False
                    }
                
                plots_by_sample[sample_id].append(plot)
                
    except FileNotFoundError:
        print(f"Error: File '{file_path}' not found.")
        sys.exit(1)
    except Exception as e:
        print(f"Error parsing plot file: {str(e)}")
        sys.exit(1)
        
    return plots_by_sample, sample_to_source

@track_performance
def add_plot_tracks(builder: CGViewBuilder, plots_by_sample: Dict[str, List[Plot]], sample_to_source: Dict[str, str], colors: List[str] = None):
    """Attach plot tracks and legend entries for each sample to the builder."""
    if colors is None:
        colors = ["rgba(0,0,255,1)", "rgba(255,0,0,1)", "rgba(0,128,0,1)", "rgba(128,0,128,1)"]
    
    color_index = 0
    for sample_id, plots in plots_by_sample.items():
        color = colors[color_index % len(colors)]
        color_index += 1
        
        for plot in plots:
            plot.meta["color"] = color
            builder.plots.append(plot)
        
        legend_item = LegendItem(
            name=sample_id,
            swatchColor=color,
            decoration="line" 
        )
        builder.legend.items.append(legend_item)
        
        source = sample_to_source[sample_id]
        
        track = Track(
            name=sample_id,
            dataType="plot",
            dataMethod="source",
            dataKeys=source,
            position="inside", 
            separateFeaturesBy="none",
            thicknessRatio=2  
        )
        builder.tracks.append(track)

@track_performance
def main():
    """Main entry, parse args, parse inputs, build a CGView JSON file."""
    parser = argparse.ArgumentParser(description='CGView Builder CLI')
    parser.add_argument('--name', default='CGViewMap', help='Name of the visualization')
    parser.add_argument('--sequence', required=True, help='Path to sequence file (FASTA format)')
    parser.add_argument('--features', help='Path to features file (GFF format)')
    parser.add_argument('--plots', help='Path to plots file (tab-delimited format with SampleID, Template, Window, Start, End, Value)')
    parser.add_argument('--output', default='CGView.json', help='Output JSON file path')
    parser.add_argument('--feature-type', default='all', help='Comma-separated list of feature types to include (default: all)')
    parser.add_argument('--feature-id', default='first', help='Attribute to use as feature ID (default: first attribute)')
    parser.add_argument('--blast', help='Path to BLAST TSV file')
    
    args = parser.parse_args()
    
    process = psutil.Process(os.getpid())
    initial_memory = process.memory_info().rss / 1024 / 1024
    print(f"Initial memory usage: {initial_memory:.2f} MB")
    
    builder = CGViewBuilder(name=args.name)
    
    if args.sequence:
        print("\nProcessing sequence file...")
        sequence = parse_sequence_file(args.sequence)
        builder.set_sequence(sequence)
    
    if args.features:
        print("\nProcessing feature file...")
        feature_types = set(args.feature_type.split(',')) if args.feature_type != 'all' else None
        features, legend_items = parse_feature_file(args.features, feature_types, args.feature_id)
        for feature in features:
            builder.add_feature(feature)

        for item in legend_items:
            builder.legend.items.append(item)
    
    if args.plots:
        print("\nProcessing plot file...")
        plots_by_sample, sample_to_source = parse_plot_file(args.plots, sequence)
        add_plot_tracks(builder, plots_by_sample, sample_to_source)
            
    if args.blast:
        print("\nProcessing BLAST file...")
        blast_features, query_to_source = parse_blast_file(args.blast)
        add_blast_tracks(builder, blast_features, query_to_source)
    
    print("\nGenerating output JSON...")
    # Serialize only after all sequence, feature, and track layers are attached to the builder.
    output = builder.to_dict()
    output_path = Path(args.output)
    
    print("\nWriting output file...")
    with output_path.open('w') as f:
        json.dump(output, f, indent=2)
        print(f"Successfully wrote output to {output_path}")
                
if __name__ == '__main__':
    main()
