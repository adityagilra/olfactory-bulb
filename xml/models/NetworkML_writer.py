#from elementtree.SimpleXMLWriter import XMLWriter # I am not using XMLWriter as it is sequential
# we want to construct a document object model dom in non-sequential order and finally write to xml
import xml.dom.ext
import xml.dom.minidom

class NetworkML_writer():
    
    def __init__(self,length_units,units):
        """
        serial order of setting sub-elements is assumed. Thus there are lots of current_<element> variables.
        However, you can simultaneously build up say instances and connections since they are sub-sub-sub-elements of disjoint populations and projections elements which are siblings and single in number.
        The populations, projections, inputs, instances elements are only created when the first population, projection, input, instance is set.
        This is because neuroml does not allow empty elements.
        units can be "SI Units" or "Physiological Units" - Note capitals.
        length_units can "meter" or "micrometer".
        """
        self.doc = xml.dom.minidom.Document()
        self.networkml = self.doc.createElement("networkml")
        self.networkml.setAttribute("xmlns","http://morphml.org/networkml/schema")
        self.networkml.setAttribute("xmlns:xsi","http://www.w3.org/2001/XMLSchema-instance")
        self.networkml.setAttribute("xmlns:meta","http://morphml.org/metadata/schema")
        self.networkml.setAttribute("xsi:schemaLocation","http://morphml.org/networkml/schema NetworkML_v1.8.1.xsd")
        self.networkml.setAttribute("length_units",length_units)
        self.units = units
        self.doc.appendChild(self.networkml)
        self.populations = None
        self.projections = None
        self.inputs = None
        
    def note(self, text):
        note_el = self.doc.createElement("meta:notes")
        self.networkml.appendChild(note_el)
        text_node = self.doc.createTextNode(text)
        note_el.appendChild(text_node)
                
    def writeToFile(self, filename):
        networkmlfile = open(filename,'w')
        xml.dom.ext.PrettyPrint(self.doc,networkmlfile)
        networkmlfile.close()
        
    # methods to specify cells:
    
    def set_population(self, name, cell_type):
        if self.populations == None:
            self.populations = self.doc.createElement("populations")
            self.networkml.appendChild(self.populations)
        self.current_population = self.doc.createElement("population")
        self.current_population.setAttribute("name",name)
        self.current_population.setAttribute("cell_type",cell_type)
        self.populations.appendChild(self.current_population)
        self.current_instances = None

    def set_instance(self, elementid,x,y,z):
        """
        size is a required attribute of instances. So everytime an instance is created, the size attribute is updated.
        """
        if self.current_instances == None:
            self.current_instances = self.doc.createElement("instances")
            self.current_instances_size = 0
            self.current_instances.setAttribute("size","0")
            self.current_population.appendChild(self.current_instances)
        instance = self.doc.createElement("instance")
        instance.setAttribute("id",elementid)
        self.current_instances.appendChild(instance)
        location = self.doc.createElement("location")
        location.setAttribute("x",x)
        location.setAttribute("y",y)
        location.setAttribute("z",z)
        instance.appendChild(location)
        self.current_instances_size += 1
        self.current_instances.setAttribute("size",str(self.current_instances_size))

    # methods to specify connections:
        
    def set_projection(self, name, source, target):
        if self.projections == None:
            self.projections = self.doc.createElement("projections")
            self.projections.setAttribute("units",self.units)
            self.networkml.appendChild(self.projections)
        self.current_projection = self.doc.createElement("projection")
        self.current_projection.setAttribute("name",name)
        self.current_projection.setAttribute("source",source)
        self.current_projection.setAttribute("target",target)
        self.projections.appendChild(self.current_projection)
        self.current_connections = None
        
    def set_synapse_props(self, synapse_type):
        synapse_props = self.doc.createElement("synapse_props")
        synapse_props.setAttribute("synapse_type",synapse_type)
        self.current_projection.appendChild(synapse_props)

    def set_connection(self, elementid,pre_cell_id,post_cell_id,weight="1"):
        """
        size is not a required attribute of connections unlike for instances.
        But we do set the size anyway. Everytime a connection is created, the size attribute is updated.
        """
        if self.current_connections == None:
            self.current_connections = self.doc.createElement("connections")
            self.current_connections_size = 0
            self.current_connections.setAttribute("size","0")
            self.current_projection.appendChild(self.current_connections)
        connection = self.doc.createElement("connection")
        connection.setAttribute("id",elementid)
        connection.setAttribute("pre_cell_id",pre_cell_id)
        connection.setAttribute("post_cell_id",post_cell_id)
        self.current_connections.appendChild(connection)
        properties = self.doc.createElement("properties")
        properties.setAttribute("weight",weight)
        connection.appendChild(properties)
        self.current_connections_size += 1
        self.current_connections.setAttribute("size",str(self.current_connections_size))
