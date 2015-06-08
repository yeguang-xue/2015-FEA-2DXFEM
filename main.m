%% Cleaning

clear
clc

nx = 11;
ny = 11;
nxel = nx -1;
nyel = ny -1;

num_nodes = nx*ny;
num_elements = nxel*nyel;

nodes_container = NodeArrayClass(num_nodes);
elements_container = ElementArrayClass(num_elements);

%%% Mesh Preprocessing

global_ls = @(X) norm(X-[5,5])-2.5;

% rho, youngs, poissons
material_1 = [1, 1, 0.3];

material_2 = [1, 10, 0.3];

% Create nodes
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
        v1_X = nodes_container.getNodeByID(conn(1)).X;
        v2_X = nodes_container.getNodeByID(conn(2)).X;
        v3_X = nodes_container.getNodeByID(conn(3)).X;
        v4_X = nodes_container.getNodeByID(conn(4)).X;
        level(1) = global_ls(v1_X);
        level(2) = global_ls(v2_X);
        level(3) = global_ls(v3_X);
        level(4) = global_ls(v4_X);
        if level>0
            element = CPS4ElementClass(nodes_container,conn,material_1);
        elseif level<0
            element = CPS4ElementClass(nodes_container,conn,material_2);
        else
            %element = CPS4ElementClass(nodes_container,conn,material_2);
            element = CPS4XElementClass(nodes_container,conn,material_1,material_2,level);
        end

        elements_container.CreateElement(id, element);
        id = id + 1;
    end
end

modelmesh = MeshClass(nodes_container,elements_container);

figure(1)
hold on;
axis equal;
xlim([0,10])
ylim([0,10])

modelmesh.plot()

theta = 0:0.01:pi*2;
xxx = 2.5*cos(theta)+5;
yyy = 2.5*sin(theta)+5;
plot(xxx,yyy,'r-','linewidth',2)

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
