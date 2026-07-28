"""
Injects the final input payload into the dedicated `scifr-data` JSON script tag inside the bundled HTML template.

Last checked by: BP
"""
import argparse
import re
import orjson

SCRIPT_TAG_PATTERN = re.compile(
    r'(<script\b[^>]*\bid="scifr-data"[^>]*\btype="application/json"[^>]*>)(.*?)(</script>)',
    re.IGNORECASE | re.DOTALL,
)


def wrap_payload(json_data):
    """Wrap the payload in the structure expected by the viewer."""
    return {"payload": json_data}


def serialize_script_payload(json_data):
    """Serialize the payload JSON for embedding safely into the HTML script tag."""
    json_bytes = orjson.dumps(json_data, option=orjson.OPT_NON_STR_KEYS)
    return json_bytes.decode("utf-8").replace("</script", "<\\/script")


def replace_scifr_data_script(content, json_data):
    """Replace the existing `scifr-data` script payload."""
    matches = list(SCRIPT_TAG_PATTERN.finditer(content))
    if not matches:
        raise ValueError("scifr-data script tag not found")
    if len(matches) > 1:
        raise ValueError("multiple scifr-data script tags found")

    match = matches[0]
    payload_text = serialize_script_payload(wrap_payload(json_data))
    return content[:match.start()] + match.group(1) + payload_text + match.group(3) + content[match.end():]


def mutate_template_memory(json_data, template_path, output_path):
    """Load the template, inject the payload, and write the final HTML report."""
    # load template dna
    with open(template_path, 'r', encoding='utf-8') as f:
        template_content = f.read()

    updated_content = replace_scifr_data_script(template_content, json_data)

    # write modified template
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(updated_content)

def mutate_report_from_file(data_json, template):
    """Load a payload from disk and mutate the template into a report file."""
    # load new json payload
    with open(data_json, 'r', encoding='utf-8') as f:
        json_data = orjson.loads(f.read().strip())

    # call the memory version
    mutate_template_memory(json_data, template, 'scifr_report.html')

def parse_arguments():
    """Parse arguments."""
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--data_json', help='path to new json payload', required=True)
    parser.add_argument('--template', help='path to template file', required=True)
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_arguments()
    mutate_report_from_file(args.data_json, args.template)
