from Bio import ExPASy, SwissProt
import xml.etree.ElementTree as ET

# Parse UniProt XML for GO terms
def extract_go_annotations(xml_file):
    go_annotations = {}
    tree = ET.parse(xml_file)
    root = tree.getroot()

    for entry in root.findall("{http://uniprot.org/uniprot}entry"):
        protein_id = entry.find("{http://uniprot.org/uniprot}accession").text
        go_terms = []

        for db_reference in entry.findall("{http://uniprot.org/uniprot}dbReference"):
            if db_reference.attrib.get("type") == "GO":
                go_id = db_reference.attrib.get("id")
                go_term = db_reference.find("{http://uniprot.org/uniprot}property[@type='term']").attrib.get("value")
                go_terms.append((go_id, go_term))

        go_annotations[protein_id] = go_terms
    return go_annotations

# Example usage
family_go_annotations = extract_go_annotations("family_sequences.xml")
swissprot_go_annotations = extract_go_annotations("swissprot.xml")

# Save or process GO annotations as needed
print(family_go_annotations)


