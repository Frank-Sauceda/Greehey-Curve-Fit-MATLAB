function [mu_hat, sigma_hat, ci99] = var_ksFit( data )

[f, xi] = ksdensity( data, 'npoints', 50 );
hist( data, xi );
hold on;
dx = xi(2)-xi(1);
plot( xi, f*length(data)*dx, 'r' );

[mf, mdx] = max( f );
mu_hat = xi( mdx );
sigma = std( data );
err_lookforward = floor(mdx + 0.5*sigma/dx);

for kk = 1:3
    sigma_hat = sigma*0.5:sigma/200:sigma*1.5;
    for i = 1:length( sigma_hat )
        y = normpdf( xi, mu_hat, sigma_hat(i) );
        my = max( y );
        delta(i) = sum((y(1:err_lookforward)./my-f(1:err_lookforward)./mf).^2);
    end
    [mx, mdx] = min( delta );
    sigma_hat = sigma_hat( mdx );
    sigma = sigma_hat;
    fprintf( 1, 'sigma: min delta = %0.3f\n', mx );
    
    mu = mu_hat;
    mu_hat = mu*0.5:mu/200:mu*1.5;
    for i = 1:length( mu_hat )
        y = normpdf( xi, mu_hat(i), sigma_hat );
        my = max( y );
        delta(i) = sum((y(1:err_lookforward)./my-f(1:err_lookforward)./mf).^2);
    end
    [mx, mdx] = min( delta );
    mu_hat = mu_hat( mdx );
    mu = mu_hat;
    fprintf( 1, 'mu: min delta = %0.3f\n', mx );
end

y = normpdf( xi, mu_hat, sigma_hat );
my = max(y);
hdl = plot( xi, normpdf( xi, mu_hat, sigma_hat ) * mf/my * length(data)*dx, 'm--' );
set( hdl, 'lineWidth', 2 );
ci99 = norminv( 0.99, mu_hat, sigma_hat );
hdl = line( [ci99 ci99], [0 mf*length(data)*dx] );
set( hdl, 'color', 'k' );
return;

