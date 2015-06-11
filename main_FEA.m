%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note: This program is just used to test standard FEA during development
% Author: Yeguang Xue (Northwestern Univ)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear, clc

%% Mesh Preprocessing
nx = 11;
ny = 11;
nxel = nx -1;
nyel = ny -1;

num_nodes = nx*ny;
num_elements = nxel*nyel;

nodes_container = NodeArrayClass(num_nodes);
elements_container = ElementArrayClass(num_elements);


rho = 1;
youngs = 1;
nu = 0.3;

material_1 = [rho, youngs, nu];

id = 1;

for i = 1:nx
    for j = 1:ny
        x = j-1;
        y = i-1;
        node = NodeClass( [x,y] );
        nodes_container.CreateNode(id, node);
        id = id + 1;
    end
end

id = 1;
for i = 1:nxel
    for j = 1:nyel
        conn = zeros(1,4);
        conn(1) = (j-1)*nx+i;
        conn(2) = (j-1)*nx+i+1;
        conn(3) = j*nx+i+1;
        conn(4) = j*nx+i;
                
        element = CPS4ElementClass(nodes_container,conn,material_1);
        elements_container.CreateElement(id, element);
        id = id + 1;
    end
end

modelmesh = MeshClass(nodes_container,elements_container);

figure(1)
hold on;
axis equal;

modelmesh.plot()

% Boundary Node Set
forceBC = [];
dispBC = zeros(1,3);
bcid = 1;
for j=1:ny
    nid = (j-1)*nx+1;
    dispBC(bcid,:)= [nid,-1,0];
    bcid = bcid +1;
    nid = j*nx;
    dispBC(bcid,:)= [nid,1,0];
    bcid = bcid +1;
end

%% Test 

K = modelmesh.build_stiff_matrix();

       
f = modelmesh.build_force_matrix(forceBC);

      
[K_modify, f_modify ] = modelmesh.applydispBC(K,f,dispBC);

d = K_modify\f_modify;

modelmesh.AssignResults(d);

modelmesh.outputVtk('test_undeformed.vtk','UNDEFORMED');
modelmesh.outputVtk('test_deformed.vtk','DEFORMED');
