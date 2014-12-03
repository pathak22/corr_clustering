function fv  = cnntest_fv(net,x)
    %  feedforward
    net = cnnff(net, x);
    fv = net.fv;
end