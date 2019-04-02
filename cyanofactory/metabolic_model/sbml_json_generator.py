"""
Copyright (c) 2018 Gabriel Kind <kind hs-mittweida de>
Hochschule Mittweida, University of Applied Sciences

Released under the MIT license
"""

import json
import metabolic_model.metabolic_model as metabolic_model

class SbmlJsonGenerator(json.JSONEncoder):
    element_filter = ["rdf:RDF", "rdf:Description", "rdf:Bag", "groups:idRef"]

    def __init__(self, *args, **kwargs):
        self.obj = {}
        self.path = []

        super().__init__(*args, **kwargs)

    def default(self, obj):
        if isinstance(obj, metabolic_model.ElementBase):
            obj.write_sbml(self)
            return self.obj
        else:
            return json.JSONEncoder.default(self, obj)

    def startPrefixMapping(self, *args, **kwargs):
        # no-op
        return

    def endPrefixMapping(self, *args, **kwargs):
        # no-op
        return

    def startElement(self, name, attrs, is_list):
        if name in SbmlJsonGenerator.element_filter:
            return

        if len(self.path) > 0:
            ref = self.obj[self.path[0]]
            for i in range(1, len(self.path)):
                if isinstance(ref, list):
                    ref = ref[-1]
                else:
                    ref = ref[self.path[i]]

        else:
            ref = self.obj

        if is_list:
            a = 0
            pass

        if name not in ref:
            if isinstance(ref, list):
                if "rdf:resource" in attrs:
                    ref.append(attrs["rdf:resource"])
                elif "groups:idRef" in attrs:
                    ref.append(attrs["groups:idRef"])
                else:
                    ref.append(attrs)
            else:
                if is_list:
                    if attrs:
                        raise ValueError("attrs not empty!")
                    ref[name] = []
                else:
                    ref[name] = attrs
        else:
            if not is_list:
                raise ValueError("Overwrite! " + name)

        self.path.append(name)

    def endElement(self, name):
        if name in SbmlJsonGenerator.element_filter:
            return

        n = self.path.pop()
        if name != n:
            raise ValueError("tag mismatch")

    def startElementNS(self, name, qname, attrs, is_list):
        if isinstance(name, list):
            name = name[1]

        if len(qname) == 0 or name == "sbml":
            fname = name
        else:
            fname = qname + ":" + name

        self.startElement(fname, attrs, is_list)

    def endElementNS(self, name, qname):
        if isinstance(name, list):
            name = name[1]

        if len(qname) == 0 or name == "sbml":
            fname = name
        else:
            fname = qname + ":" + name

        self.endElement(fname)

    class ElementScopeHandler(object):
        def __init__(self, parent, name, attrs=None, is_list=False):
            if attrs is None:
                attrs = {}

            self.parent = parent
            self.name = name
            self.attrs = attrs
            self.is_list = is_list

        def make_attr(self):
            attr_out = []

            for i in range(0, len(self.attrs), 2):
                val = self.attrs[i + 1]
                if len(str(val)) == 0:
                    continue
                attr_out.append([self.attrs[i], val])

            return dict(attr_out)

        def __enter__(self):
            self.parent.startElement(self.name, self.make_attr(), self.is_list)

        def __exit__(self, exception_type, exception_value, traceback):
            self.parent.endElement(self.name)

    class ElementScopeHandlerNS(object):
        def __init__(self, parent, name, qname, attrs=None, is_list=False):
            if attrs is None:
                attrs = {}

            self.parent = parent
            self.name = name
            self.attrs = attrs
            self.qname = qname
            self.is_list = is_list

        def make_attr_ns(self):
            attr_out = []

            for i in range(0, len(self.attrs), 3):
                ns, name, val = self.attrs[i:i + 3]

                if len(str(val)) == 0:
                    continue

                if len(ns) == 0:
                    attr_out.append([name, val])
                else:
                    attr_out.append([ns + ":" + name, val])

            return dict(attr_out)

        def __enter__(self):
            self.parent.startElementNS(self.name, self.qname, self.make_attr_ns(), self.is_list)

        def __exit__(self, exception_type, exception_value, traceback):
            self.parent.endElementNS(self.name, self.qname)

    def element(self, name, attrs=None, is_list=False):
        return SbmlJsonGenerator.ElementScopeHandler(self, name, attrs, is_list)

    def elementNS(self, ns, name, attrs=None, is_list=False):
        return SbmlJsonGenerator.ElementScopeHandlerNS(self, name, ns, attrs, is_list)
