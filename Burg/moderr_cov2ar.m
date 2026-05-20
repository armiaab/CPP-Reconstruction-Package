function model_error =  moderr_cov2ar(cov,truear,truema,N,long_order)
%function model_error =  moderr_long_ar(ar,ma,truear,truema,N,long_order)
% Computes model error by transforming the covariance function cov to an
% AR(long_order) process with levinson.
cov=cov(:)';
if length(cov)<long_order+1
    cov=[cov zeros(1,long_order+1-length(cov))];
end
par=levinson(cov,long_order);
model_error=moderr(par,1,truear,truema,N);
end

