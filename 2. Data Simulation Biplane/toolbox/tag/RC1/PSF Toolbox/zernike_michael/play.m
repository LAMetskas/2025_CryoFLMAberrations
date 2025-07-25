clear all
close all

Z = Zernike_Polynomials;
Z.setN_Wyant(5);
Z.initialize_Wyant();
Z.setN_Noll(5);
Z.initialize_Noll();
Z

Z.displayC_Wyant();
Z.displayC_Noll();

fprintf('\n');
l = Z.indices_WyantNMtoL(3, 2);
fprintf('l = %d\n', l);
[n, m] = Z.indices_WyantLtoNM(11);
fprintf('(n, m) = (%d, %d)\n', n, m);
l = Z.indices_WyantNMtoL(3, -2);
fprintf('l = %d\n', l);
[n, m] = Z.indices_WyantLtoNM(12);
fprintf('(n, m) = (%d, %d)\n', n, m);
l = Z.indices_WyantNMtoL(3, 0);
fprintf('l = %d\n', l);
[n, m] = Z.indices_WyantLtoNM(15);
fprintf('(n, m) = (%d, %d)\n', n, m);

fprintf('\n');
l = Z.indices_NollNMtoL(3, 3);
fprintf('l = %d\n', l);
[n, m] = Z.indices_NollLtoNM(10);
fprintf('(n, m) = (%d, %d)\n', n, m);
l = Z.indices_NollNMtoL(3, -3);
fprintf('l = %d\n', l);
[n, m] = Z.indices_NollLtoNM(9);
fprintf('(n, m) = (%d, %d)\n', n, m);
l = Z.indices_NollNMtoL(4, 0);
fprintf('l = %d\n', l);
[n, m] = Z.indices_NollLtoNM(11);
fprintf('(n, m) = (%d, %d)\n', n, m);

fprintf('\n');
fprintf('Wyant: Z_3^1(1, 1) = %f\n', Z.poly_WyantNM(3, 1, 1, 1));
fprintf('Wyant: Z_3^1(2, 2) = %f\n', Z.poly_WyantNM(3, 1, 2, 2));
fprintf('Wyant: Z_3^1(3, 3) = %f\n', Z.poly_WyantNM(3, 1, 3, 3));
fprintf('Wyant: Z_3^1(4, 4) = %f\n', Z.poly_WyantNM(3, 1, 4, 4));
Z.poly_WyantNM(3, 1, [1, 2; 3, 4], [1, 2; 3, 4])

fprintf('Noll:  Z_3^1(1, 1) = %f\n', Z.poly_NollNM(3, 1, 1, 1));
fprintf('Noll:  Z_3^1(2, 2) = %f\n', Z.poly_NollNM(3, 1, 2, 2));
fprintf('Noll:  Z_3^1(3, 3) = %f\n', Z.poly_NollNM(3, 1, 3, 3));
fprintf('Noll:  Z_3^1(4, 4) = %f\n', Z.poly_NollNM(3, 1, 4, 4));
Z.poly_NollNM(3, 1, [1, 2; 3, 4], [1, 2; 3, 4])

fprintf('Wyant: Z_0^0 (1, 1) = %f\n', Z.poly_WyantNM(0,  0, 1, 1));
fprintf('Wyant: Z_1^0 (1, 1) = %f\n', Z.poly_WyantNM(1,  0, 1, 1));
fprintf('Wyant: Z_1^1 (1, 1) = %f\n', Z.poly_WyantNM(1,  1, 1, 1));
fprintf('Wyant: Z_1^-1(1, 1) = %f\n', Z.poly_WyantNM(1, -1, 1, 1));

fprintf('Noll:  Z_0^1 (1, 1) = %f\n', Z.poly_NollNM(0,  0, 1, 1));
fprintf('Noll:  Z_1^1 (1, 1) = %f\n', Z.poly_NollNM(1,  1, 1, 1));
fprintf('Noll:  Z_1^-1(1, 1) = %f\n', Z.poly_NollNM(1, -1, 1, 1));
fprintf('Noll:  Z_2^0 (1, 1) = %f\n', Z.poly_NollNM(2,  0, 1, 1));

Z.setN_Wyant(3);
nZ_terms = Z.n_Zernikes_Wyant(Z.Wnmax); 

fprintf('Wyant: %f\n', Z.poly_sum_Wyant(1, 1, 1));
fprintf('Wyant: %f\n', Z.poly_sum_Wyant(2, 2, 1));
fprintf('Wyant: %f\n', Z.poly_sum_Wyant(3, 3, 1));
fprintf('Wyant: %f\n', Z.poly_sum_Wyant(4, 4, 1));
Z.poly_sum_Wyant([1, 2; 3, 4], [1, 2; 3, 4], ones(2, 2), ones(1, nZ_terms))
Z.poly_sum_Wyant([1, 2; 3, 4], [1, 2; 3, 4], ones(2, 2))
Z.poly_sum_Wyant([1, 2; 3, 4], [1, 2; 3, 4])

Z.setN_Noll(3);
nZ_terms = Z.n_Zernikes_Noll(Z.Wnmax); 

fprintf('Noll:  %f\n', Z.poly_sum_Noll(1, 1, 1));
fprintf('Noll:  %f\n', Z.poly_sum_Noll(2, 2, 1));
fprintf('Noll:  %f\n', Z.poly_sum_Noll(3, 3, 1));
fprintf('Noll:  %f\n', Z.poly_sum_Noll(4, 4, 1));
Z.poly_sum_Noll([1, 2; 3, 4], [1, 2; 3, 4], ones(2, 2), ones(1, nZ_terms))
Z.poly_sum_Noll([1, 2; 3, 4], [1, 2; 3, 4], ones(2, 2))
Z.poly_sum_Noll([1, 2; 3, 4], [1, 2; 3, 4])

% -----------------------------------------------------------------------------

Z = Zernike_Polynomials;

for i = 1 : 2

   if i == 1
      Z.Ordering = 'Wyant';
   else
      Z.Ordering = 'Noll';
   end

   Z.setN(5);
   Z.initialize();
   Z

   Z.displayC();
   fprintf('\n');

   Z.setN(3);
   nZ_terms = Z.n_Zernikes(3); 

   fprintf('%s: %f\n', Z.Ordering, Z.poly_sum(1, 1, 1));
   fprintf('%s: %f\n', Z.Ordering, Z.poly_sum(2, 2, 1));
   fprintf('%s: %f\n', Z.Ordering, Z.poly_sum(3, 3, 1));
   fprintf('%s: %f\n', Z.Ordering, Z.poly_sum(4, 4, 1));
   Z.poly_sum([1, 2; 3, 4], [1, 2; 3, 4], ones(2, 2), ones(1, nZ_terms))
   Z.poly_sum([1, 2; 3, 4], [1, 2; 3, 4], ones(2, 2))
   Z.poly_sum([1, 2; 3, 4], [1, 2; 3, 4])

   ZM = Z.matrix_Z([1, 2; 3, 4], [1, 2; 3, 4]);
   size(ZM)
   sum(ZM, 3)

end
