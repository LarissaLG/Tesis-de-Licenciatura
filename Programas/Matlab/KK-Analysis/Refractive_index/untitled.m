% Define the file path
filePath = 'C:\Users\Poh\Documents\GitHub\Tesis-de-Licenciatura\Programas\ParteIm.csv';

% Read the data into a table, skipping the first 3 rows
dataTable = readtable(filePath, 'HeaderLines', 1);

%%
omega = dataTable{1:112, 1};         % eV
k      = dataTable{1:112, 2}; % segunda columna, distinto rango

%%
%Aplicar KK
nKK = kkrebook(omega, k);


subplot(2,1,1);
x = omega;
y1 = nKK;
plot(x,y1, "r--")

subplot(2,1,2); 
y2 = k;
plot(x,y2)
