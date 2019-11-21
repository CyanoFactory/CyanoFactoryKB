"""
Copyright (c) 2018 Gabriel Kind <kind hs-mittweida de>
Hochschule Mittweida, University of Applied Sciences

Released under the MIT license
"""

from xml.sax.saxutils import XMLGenerator
from xml.sax.xmlreader import AttributesNSImpl

class SbmlXMLGenerator(XMLGenerator):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    class ElementScopeHandler(object):
        def __init__(self, parent, name, attrs=None):
            if attrs is None:
                attrs = {}

            self.parent = parent
            self.name = name
            self.attrs = attrs

        def make_attr(self):
            attr_out = []

            for i in range(0, len(self.attrs), 2):
                val = self.attrs[i + 1]
                if isinstance(val, bool):
                    val = str(val).lower()
                else:
                    val = str(val)

                if len(val) == 0:
                    continue

                attr_out.append([self.attrs[i], val])

            return dict(attr_out)

        def __enter__(self):
            self.parent.startElement(self.name, self.make_attr())

        def __exit__(self, exception_type, exception_value, traceback):
            self.parent.endElement(self.name)

    class ElementScopeHandlerNS(object):
        namespaces = dict(
            fbc="http://www.sbml.org/sbml/level3/version1/fbc/version2",
            groups="http://www.sbml.org/sbml/level3/version1/groups/version1",
            rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#",
            bqbiol="http://biomodels.net/biology-qualifiers/",
            bqmodel="http://biomodels.net/model-qualifiers/",
            wedesign="https://cyanofactory.hs-mittweida.de/wedesign/version1"
        )

        def __init__(self, parent, name, qname, attrs=None):
            if attrs is None:
                attrs = {}

            if len(qname) == 0:
                name = ["", name]
            elif qname in SbmlXMLGenerator.ElementScopeHandlerNS.namespaces:
                name = [SbmlXMLGenerator.ElementScopeHandlerNS.namespaces[qname], name]

            self.parent = parent
            self.name = name
            self.attrs = attrs
            self.qname = qname

        def make_attr_ns(self):
            arg1 = {}
            arg2 = {}

            for i in range(0, len(self.attrs), 3):
                ns, name, val = self.attrs[i:i+3]

                if isinstance(val, bool):
                    val = str(val).lower()
                else:
                    val = str(val)

                if len(val) == 0:
                    continue

                if len(ns) == 0:
                    ns_uri = ""
                else:
                    ns_uri = SbmlXMLGenerator.ElementScopeHandlerNS.namespaces[ns]

                arg1[(ns_uri, name)] = val
                arg2[(ns_uri, name)] = ns + ":" + name

            return AttributesNSImpl(arg1, arg2)

        def __enter__(self):
            self.parent.startElementNS(self.name, self.qname, self.make_attr_ns())

        def __exit__(self, exception_type, exception_value, traceback):
            self.parent.endElementNS(self.name, self.qname)

    def element(self, name, attrs=None, is_list=False):
        return SbmlXMLGenerator.ElementScopeHandler(self, name, attrs)

    def elementNS(self, ns, name, attrs=None, is_list=False):
        return SbmlXMLGenerator.ElementScopeHandlerNS(self, name, ns, attrs)

class SbmlXMLGeneratorWithWeDesign(SbmlXMLGenerator):
    def __init__(self, *args, **kwargs):
        self.enable_wedesign = True

        super().__init__(*args, **kwargs)
