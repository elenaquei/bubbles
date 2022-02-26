function c_cell = compose_cell(A11, A12, A13, A21, A22, A23, A31, A32, A33)
% INPUT
% 9 cells of same size
% OUTPUT
% bigger cell including all the subcells

SizeA = size(A11);
full_size = 3*SizeA;
c_cell = cell(full_size);
for i = 1:SizeA(1)
    for j = 1:SizeA(2)
        c_cell{i,j} = A11{i,j};
        c_cell{i,SizeA(2)+j} = A12{i,j};
        c_cell{i,SizeA(2)*2+j} = A13{i,j};
        c_cell{SizeA(1)+i,j} = A21{i,j};
        c_cell{SizeA(1)+i,SizeA(2)+j} = A22{i,j};
        c_cell{SizeA(1)+i,SizeA(2)*2+j} = A23{i,j};
        c_cell{2*SizeA(1)+i,j} = A31{i,j};
        c_cell{2*SizeA(1)+i,SizeA(2)+j} = A32{i,j};
        c_cell{2*SizeA(1)+i,SizeA(2)*2+j} = A33{i,j};
    end
end