function Y = zgPCATrans( X, Dic, Xm )

Xc = X-repmat(Xm,1,size(X,2)); 

Y = Dic.'*Xc;

end
