function Ycont = YHopfCont(A0, A1,X0,X1,G0,G1,const_G,Es0,Es1, const_Es0, const_Es1)

A0_op = compose_overoperator(A0,X0.size_real, 2*X0.node_space,2*X0.node_time);
A1_op = compose_overoperator(A1,X1.size_real, 2*X1.node_space,2*X1.node_time);

    function Y_bound=max_YboundHopf()
        % Y at the extrema
        Y0=...
            A0_op*...
            evaluate_HopfKS(X0,G0,const_G,Es0,const_Es0,1);
        Y1=...
            A1_op*...
            evaluate_HopfKS(X1,G1,const_G,Es1,const_Es1,1);
        
        Y_bound = norm(max(abs(Y0),abs(Y1)));
        
    end

    function CB = centerYHopf()
        
        ADelta = A0-A1;
        %ADelta_op = A0_op - A1_op;
        x0 = PDE_vec2vec(X0);
        x1 = PDE_vec2vec(X1);
        xDelta = x0-x1;
        
        xs_real = infsup(min(real(x0),real(x1)),max(real(x0),real(x1)));
        xs_imag = infsup(min(imag(x0),imag(x1)),max(imag(x0),imag(x1)));
        xs = xs_real + 1i* xs_imag;
        
        Xs = PDE_vector(xs, 3, 1, X0.node_space);
        XDelta = PDE_vector(xDelta, 3, 1, X0.node_space);
        cDelta = const_Es0-const_Es1;
        
        DDH = 0*intval(PDE_vec2vec(reshape(X0,X0.node_space*2,X0.node_time*2)));
        DDH(1) = xs * xDelta.'+cDelta;
        DDH_vec = 0*intval(X0);
        DDH_vec.parameters(1) =  xs * xDelta.'+cDelta;
        first_term = norm(PDE_vector(2*ADelta*DDH.', 3, 1, 2*X0.node_space));
        %first_term_2 = norm(2*ADelta_op*DDH_vec);
        
        Es_real = infsup(min(real(Es0),real(Es1)),max(real(Es0),real(Es1)));
        Es_imag = infsup(min(imag(Es0),imag(Es1)),max(imag(Es0),imag(Es1)));
        Es_s = Es_real + 1i* Es_imag;
        
        Gs_real = infsup(min(real(G0),real(G1)),max(real(G0),real(G1)));
        Gs_imag = infsup(min(imag(G0),imag(G1)),max(imag(G0),imag(G1)));
        G_s = Gs_real + 1i* Gs_imag;
        %Es_s = infsup(min(Es0,Es1),max(Es0,Es1));
        DxHxs = derivative_Hopf(Xs, G_s, Es_s, 1);
        
        DxHxs_xDelta_long = DxHxs*PDE_vec2vec(reshape(XDelta, X0.node_space*2, X0.node_time*2)).';
        
        second_term = norm(PDE_vector(ADelta*DxHxs_xDelta_long, 3, 1, X0.node_space*2));
        
        As = infsup(min(real(A0),real(A1))+1i*min(imag(A0),imag(A1)), ...
            max(real(A0),real(A1))+1i*max(imag(A0),imag(A1)));
       % DDDHxDxD = second_der(Xs, XDelta);
        
        DDDHxDxD = (second_derHopf(Xs, XDelta,1)*...
            PDE_vec2vec(reshape(XDelta, XDelta.node_space*2, XDelta.node_time*2)).');
        
        third_term = norm(PDE_vector(As*DDDHxDxD, 3, 1, X0.node_space*2));
        
        CB = first_term+second_term+third_term;
    end

max_bounds =max_YboundHopf();
center_bound = centerYHopf();
Ycont = max_bounds + 1/8* center_bound;
end





% function DDDHxDxD = second_der(Xs, V)
% % INPUT
% % Xs    PDE_vector
% % V     PDE_vector
% % OUTPUT
% % DDDHxDxD  PDE_vector
% 
% Xs = reshape(Xs, Xs.node_space*2, Xs.node_time*2);
% V = reshape(V, V.node_space*2, V.node_time*2);
% 
% DDDHxDxD = Xs;
% DDDHxDxD.parameters = 0*DDDHxDxD.parameters;
% 
% K2_y = derivative_operator(2,Xs.vector1D);
% K2_z = derivative_operator(2,squeeze(Xs.vector2D));
% 
% ys = Xs.vector1D;
% zs = squeeze(Xs.vector2D);
% l1s = Xs.parameters(1);
% l2s = Xs.parameters(2);
% as = Xs.parameters(3);
% 
% y = V.vector1D;
% z = squeeze(V.vector2D);
% l1 = V.parameters(1);
% l2 = V.parameters(2);
% a = V.parameters(3);
% 
% DDDHxDxD.vector1D(:) = -K2_y.^4*(2*ys*l1*l2+2*l2s*l1*y+2*y*l2*l1).'+2*K2_y.^2*y.'*l1+...
%     4*1i*K2_y*(2*l1*Convo(y,ys)+l1*Convo(y,y)).';
% 
% 
% zs = reshape(zs,1,[]);
% z = reshape(z,1,[]);
% y = reshape(y,1,[]);
% 
% DDDHxDxD.vector2D(:,:) = -K2_z.^4*(2*zs*l1*l2+2*l2s*l1*z+2*z*l2*l1s).'+K2_z.^4*2*z.'*l1+...
%     4*1i*K2_z*(2*l1*Convo(zs,y)+2*l1s*Convo(z,y)+2*l1*Convo(z,ys)).'+...
%     4*1i*K2_z*(l1*(a*2*Convo(zs,zs)+4*as*Convo(zs,z))+l1s*(4*a*Convo(z,zs)+2*as*Convo(z,z))).';
% 
% end





