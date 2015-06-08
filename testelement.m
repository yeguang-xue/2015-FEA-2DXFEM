clc,clear

nodes_container = NodeArrayClass(4);

node = NodeClass( [0,1] );
nodes_container.CreateNode(1, node);

node = NodeClass( [0,0] );
nodes_container.CreateNode(2, node);

node = NodeClass( [2,0.5] );
nodes_container.CreateNode(3, node);

node = NodeClass( [2,1] );
nodes_container.CreateNode(4, node);

rho = 1;
youngs = 3E7;
nu = 0.3;

material_1 = [rho, youngs, nu];
conn = [1,2,3,4];

element = CPS4ElementClass(nodes_container,conn,material_1);

element.derative_matrix(element.gauss(1,:));
element.build_stiff_matrix();


elements_container = ElementArrayClass(1);
elements_container.CreateElement(1, element);


modelmesh = MeshClass(nodes_container,elements_container);
K = modelmesh.build_stiff_matrix();

forceBC = [1, 0, -20;
           4, 0, -20];
       
f = modelmesh.build_force_matrix(forceBC);

dispBC = [1, 0, 0;
          2, 0, 0];
      
[K_modify, f_modify ] = modelmesh.applydispBC(K,f,dispBC)

K_modify\f_modify