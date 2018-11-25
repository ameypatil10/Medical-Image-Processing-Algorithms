function ans = myFilter(Rf, theta, filter_name, rect_width)
    rect = @(w,L) [(abs(w) <= L)];
    filter_names = {'Ram-Lak', 'Shepp-Logan', 'Cosine'};
    filter_definitions = {@(w,L) abs(w) .* rect(w,L),
        @(w,L) abs(w) .* ((w ~= 0).*sin((pi.*w) /(2*L)) + (w==0))  ./ ((w ~= 0).*(pi.*w) /(2*L)+ (w==0)) .* rect(w,L),
        @(w,L) abs(w) .* cos(pi*w/(2*L)) .* rect(w,L)};
    function_map = containers.Map(filter_names,filter_definitions);
    if ~any(strcmp(filter_names,filter_name))
        A = @(w,L) 1;
    else
        A = function_map(filter_name);
    end
    freqs = linspace(-1, 1, size(Rf,1)).';
    FRf = fftshift(fft(Rf),1);
    wFRf = FRf .* repmat(A(freqs,rect_width),[1 size(Rf,2)]);
    wFRf = ifftshift(wFRf,1);
    F1wFRf = real(ifft(wFRf,[],1));
    ans = iradon(F1wFRf, theta, 'nearest', 'none', 1.0, 256);
end