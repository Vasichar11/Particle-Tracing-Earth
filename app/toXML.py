import xml.etree.ElementTree as ET
import sys

def convert_header_to_xml(file_path):
    with open(file_path, 'r') as file:
        header_content = file.read()

    constants = ET.Element('constants')

    lines = header_content.split('\n')
    for line in lines:
        if line.startswith("const"):
            parts = line.split("=")
            name = parts[0].split()[-1]
            value = parts[1].split()[0]
            description = line.split("//")[-1].strip()

            parameter = ET.SubElement(constants, 'parameter')
            name_elem = ET.SubElement(parameter, 'name')
            value_elem = ET.SubElement(parameter, 'value')
            description_elem = ET.SubElement(parameter, 'description')

            name_elem.text = name
            value_elem.text = value
            description_elem.text = description

    xml_content = ET.tostring(constants, encoding='utf-8', method='xml').decode()

    with open('constants.xml', 'w') as output_file:
        output_file.write(xml_content)


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Please provide the header file path as a command line argument.")
        print("Usage: python script.py <header_file>")
        sys.exit(1)

    header_file = sys.argv[1]
    convert_header_to_xml(header_file)
