from .sbml_parser import *
import xml.dom.minidom

# Parse BIGG SBML
sbml_handler = SbmlHandler()
push_handler(sbml_handler)
parser.parse("iJN678.xml")

def sbml_read_write():
    # TODO
    #OptGeneParser.from_model(sbml_handler.model).write(open("Toy_model_VLCsubject.txt", "w"))

    # Write SBML
    xml_handle = StringIO()
    writer = SbmlXMLGenerator(xml_handle, "utf-8")

    writer.startDocument()
    sbml_handler.model.write_sbml(writer)
    writer.endDocument()

    dom = xml.dom.minidom.parseString(xml_handle.getvalue())
    out = open("iJN678_out.xml", "w")
    out.write(dom.toprettyxml())

    # Parse our SBML & write it
    sbml_handler = SbmlHandler()
    push_handler(sbml_handler)
    parser.parse("iJN678_out.xml")

    xml_handle = StringIO()
    writer = SbmlXMLGenerator(xml_handle, "utf-8")

    writer.startDocument()
    sbml_handler.model.write_sbml(writer)
    writer.endDocument()

    dom = xml.dom.minidom.parseString(xml_handle.getvalue())
    out = open("iJN678_out2.xml", "w")
    out.write(dom.toprettyxml())

    # Write JSON
    out = open("iJN678_out.json", "w")
    json.dump(sbml_handler.model, out, cls=SbmlJsonGenerator, indent='\t')
    js = json.dumps(sbml_handler.model, cls=SbmlJsonGenerator)

    # Read & Write our JSON
    mm = MetabolicModel.from_json(json.loads(js)["sbml"]["model"])
    out = open("iJN678_out2.json", "w")
    json.dump(mm, out, cls=SbmlJsonGenerator, indent='\t')

    a = 0

