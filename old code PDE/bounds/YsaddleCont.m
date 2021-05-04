function Ycont = YsaddleCont(A0, A1,X0,X1,G,const_G,Es0,Es1, const_Es0, const_Es1,EDelta, const_EDelta)

A0_op = compose_overoperator3(A0,X0.size_real, 2*X0.node_space,2*X0.node_time);
A1_op = compose_overoperator3(A1,X1.size_real, 2*X1.node_space,2*X1.node_time);

    function Y_bound=max_Yboundsaddle()
        % Y at the extrema
        y0 = evaluate_saddle(X0,G,const_G,Es0,const_Es0,EDelta, const_EDelta,1);
        Y0= A0_op*y0;
        
        y1 = evaluate_saddle(X1,G,const_G,Es1,const_Es1,EDelta, const_EDelta,1);
        Y1= A1_op*y1;
        
        Y_bound = norm(max(abs(Y0),abs(Y1)));
        
    end

    function CB = centerYsaddle()
        
        ADelta = A0-A1;
        %ADelta_op = A0_op - A1_op;
        x0 = PDE_vec2vec(X0);
        x1 = PDE_vec2vec(X1);
        xDelta = x0-x1;
        
        xs_real = infsup(min(real(x0),real(x1)),max(real(x0),real(x1)));
        xs_imag = infsup(min(imag(x0),imag(x1)),max(imag(x0),imag(x1)));
        xs = xs_real + 1i* xs_imag;
        
        Xs = PDE_vector(xs, 9, 3, X0.node_space);
        XDelta = PDE_vector(xDelta, 9, 3, X0.node_space);
        cDelta = const_Es0-const_Es1;
        
        DDH = 0*intval(PDE_vec2vec(reshape(X0,X0.node_space*2,X0.node_time*2)));
        DDH(1) = xs * xDelta.'+cDelta;
        DDH_vec = 0*intval(X0);
        DDH_vec.parameters(1) =  xs * xDelta.'+cDelta;
        Y = PDE_vector(2*ADelta*DDH.', 9, 3, 2*X0.node_space);
        first_term = norm(PDE_vector(2*ADelta*DDH.', 9, 3, 2*X0.node_space));
        %first_term_2 = norm(2*ADelta_op*DDH_vec);
        
        Es_real = infsup(min(real(Es0),real(Es1)),max(real(Es0),real(Es1)));
        Es_imag = infsup(min(imag(Es0),imag(Es1)),max(imag(Es0),imag(Es1)));
        Es_s = Es_real + 1i* Es_imag;
        
        DxHxs = derivative_saddle(Xs, G, Es_s,EDelta, 1);
        
        DxHxs_xDelta_long = DxHxs*PDE_vec2vec(reshape(XDelta,X0.node_space*2,X0.node_time*2)).';
        
        second_term = norm(PDE_vector(ADelta*DxHxs_xDelta_long.', 9, 3, X0.node_space*2));
        
        As = infsup(min(real(A0),real(A1))+1i*min(imag(A0),imag(A1)), ...
            max(real(A0),real(A1))+1i*max(imag(A0),imag(A1)));
       % DDDHxDxD = second_der(Xs, XDelta);
        
       DDHxD = second_der_saddle(Xs, XDelta,1);
        DDHxDxD = PDE_vec2vec(DDHxD*...
            reshape(XDelta, XDelta.node_space*2, XDelta.node_time*2)).';
        
        third_term = norm(PDE_vector(As*DDHxDxD, 9, 3, X0.node_space*2));
        
        CB = first_term+second_term+third_term;
    end

max_bounds =max_Yboundsaddle();
center_bound = centerYsaddle();
Ycont = max_bounds + 1/8* center_bound;
end

