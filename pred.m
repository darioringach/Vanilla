function pred = pred(x)

global X ds

ncell = size(X,2);
pred = zeros(size(X));

% generate prediction for each cell

for col=1:ncell
    if(ds~=5)
        pred(:,col) = pred_one(x,col);
    else
        pred(:,col) = pred_oned(x,col);
    end
end


function y = pred_one(x,col)

global X

ncell = size(X,2);

s = x(1);   % sigma
th = x(2);  % theta
b = x(3);   % beta
a = x(4);   % alpha

% Deal with NaNs at the end of the record

y = X(:,col);
nsamp = find(isnan(y),1)-1;
if(isempty(nsamp))
    nsamp = size(X,1);
end

% Odd filter

t = 0:nsamp-1;
w = t.*exp(-t.^2 / (2*s^2));
w(2:end) = (w(2:end)-w(end:-1:2));
w = -w';
w = w/norm(w);

% Even filter

w0 = zeros(nsamp,1);
w0 = exp(-t.^2 / (2*s^2));
w0(2:end) = (w0(2:end)+w0(end:-1:2));
w0 = w0';
w0 = w0/norm(w0);
 
% Filtered signals

wf0 = fft(w0);
xf0 = real(ifft(fft(y(1:nsamp)).*wf0));
xf0 = zscore(xf0);

wf = fft(w);
xf = real(ifft(fft(y(1:nsamp)).*wf));
xf = zscore(xf);

% Of course one can combine the filters first and convolve once...
% but for historical reasons I kept them separate.

% Linear combination of filtered signals

xf = cosd(a)*xf+sind(a)*xf0;

% Output nonlinearity

y(1:nsamp) = (xf-th).^b .* (xf>=th);


% Same as above but allows for an extra delay parameter...

function y = pred_oned(x,col)

global X

ncell = size(X,2);

s = x(1);   
th = x(2);
b = x(3);
a = x(4);
d = x(5);   % delay

y = X(:,col);
nsamp = find(isnan(y),1)-1;
if(isempty(nsamp))
    nsamp = size(X,1);
end

t = -50:50;
h = (t-d).*exp(-(t-d).^2 / (2*s^2));
h = -h;
h = h/norm(h);

w = zeros(1,nsamp);
w(1:51) = h(51:end);
w(end-49:end) = h(1:50);


h = exp(-(t-d).^2 / (2*s^2));
h = h/norm(h);

w0 = zeros(1,nsamp);
w0(1:51) = h(51:end);
w0(end-49:end) = h(1:50);

w0 = w0';
w = w';

w0 = w0/norm(w0);
w = w/norm(w);
 
w = w - (w' * w0) * w0;
w = w/norm(w);

wf0 = fft(w0);
xf0 = real(ifft(fft(y(1:nsamp)).*wf0));
xf0 = zscore(xf0);

wf = fft(w);
xf = real(ifft(fft(y(1:nsamp)).*wf));
xf = zscore(xf);

xf = cosd(a)*xf+sind(a)*xf0;

y(1:nsamp) = (xf-th).^b .* (xf>=th);

