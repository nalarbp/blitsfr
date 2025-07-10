"""
Mutator of SCIFR template
"""
import argparse
import orjson

def find_and_replace_json_block(content, startIdx, endIdx, new_json_data):
    start_marker = f'JSON.parse(\'{{"startIdx":"{startIdx}"'
    end_marker = f'"endIdx":"{endIdx}"}}\')' 
    
    #find start
    start_pos = content.find(start_marker)
    if start_pos == -1:
        raise ValueError(f"startIdx not found: startIdx '{startIdx}'")
    
    #find end
    end_pos = content.find(end_marker, start_pos)
    if end_pos == -1:
        raise ValueError(f"endIdx not found: endIdx '{endIdx}'")
    
    #calculate full end position including the marker
    full_end_pos = end_pos + len(end_marker)
    
    #prepare new json payload with proper escaping
    json_bytes = orjson.dumps(new_json_data, option=orjson.OPT_NON_STR_KEYS)
    new_json_str = json_bytes.decode('utf-8').replace('\\', '\\\\').replace("'", "\\'")
    
    #construct replacement block
    replacement = f"JSON.parse('{new_json_str}')"
    
    #perform molecular cutting and ligation
    new_content = content[:start_pos] + replacement + content[full_end_pos:]
    
    return new_content

def mutate_template_memory(json_data, template_path, output_path, startIdx, endIdx):
    #load template dna
    with open(template_path, 'r', encoding='utf-8') as f:
        template_content = f.read()
    
    #perform molecular replacement
    try:
        updated_content = find_and_replace_json_block(template_content, startIdx, endIdx, json_data)
        print(f"successful cutting at idx markers: {startIdx} -> {endIdx}")
    except ValueError as e:
        raise ValueError(f"idx recognition failed: {e}")
    
    #write modified template
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(updated_content)

def mutate_report_from_file(data_json, template, startIdx, endIdx):
    #load new json payload
    with open(data_json, 'r', encoding='utf-8') as f:
        json_data = orjson.loads(f.read().strip())
    
    #call the memory version
    mutate_template_memory(json_data, template, 'scifr_report.html', startIdx, endIdx)

# CLI interface remains the same
def parse_arguments():
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--data_json', help='path to new json payload', required=True)
    parser.add_argument('--template', help='path to template file', required=True)
    parser.add_argument('--startIdx', help='tag for startIdx recognition', required=True)
    parser.add_argument('--endIdx', help='tag for endIdx recognition', required=True)
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_arguments()
    mutate_report_from_file(args.data_json, args.template, args.startIdx, args.endIdx)
    print("Mutation complete!")