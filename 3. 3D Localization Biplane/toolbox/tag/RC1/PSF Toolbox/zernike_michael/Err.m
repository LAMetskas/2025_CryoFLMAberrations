function e = Err(data, model)

   %[nx, ny, nz] = size(data);
   e = sum(abs(data(:) - model(:)))
   e = sum((data(:) - model(:)).^2 ./ model(:))
   e = mse(data, model, [])

end
