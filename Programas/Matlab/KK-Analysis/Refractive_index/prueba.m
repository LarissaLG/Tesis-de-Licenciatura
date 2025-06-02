function rechi=prueba(omega,imchi,alpha)
    %The program inputs are 1) omega, vector of the frequency
    %(or energy) components, 2) imchi, vector of the imaginary
    %part of the susceptibility under examination, and 3) alpha,
    %the value of the moment considered. The two vectors
    %1) and 2) must have the same length.
    %The output is the estimate of the real part as obtained
    %with K-K relations.
    %In order to use this program, save the whole text contained
    %in this section in a file and name it kkrebook.m

    if size(omega,1)>size(omega,2)
        omega=omega';
    end
    if size(imchi,1)>size(imchi,2)
    imchi=imchi';
    end
    %Here the program rearranges the two vectors so that,
    %whichever their initial shape, they become row vectors.
    g=size(omega,2);
    %Size of the vectors.%
    rechi=zeros(size(imchi));
    %The output is initialized.

    a=zeros(size(imchi));
    b=zeros(size(imchi));
    %Two vectors for intermediate calculations are initialized

    deltaomega=omega(2)-omega(1);
    %Here we compute the frequency (or energy) interval
    j=1;
    beta1=0;
    for k=2:g
        b(1)=beta1+imchi(k)*omega(k)^(2*alpha+1)/(omega(k)^2-omega(1)^2);
        beta1=b(1);
    end

    rechi(1)=2/pi*deltaomega*b(1)*omega(1)^(-2*alpha);
    %First element of the output: the principal part integration
    %is computed by excluding the first element of the input
    j=g;
    alpha1=0;
    for k=1:g-1
        a(g)=alpha1+imchi(k)*omega(k)^(2*alpha+1)/(omega(k)^2-omega(g)^2);
        alpha1=a(g);
    end

    rechi(g)=2/pi*deltaomega*a(g)*omega(g)^(-2*alpha);
    %Last element of the output: the principal part integration
    %is computed by excluding the last element of the input
    
    for j=2:g-1
    %Loop on the inner components of the output vector.
        alpha1=0;
        beta1=0;
        for k=1:j-1
            a(j)=alpha1+imchi(k)*omega(k)^(2*alpha+1)/(omega(k)^2-omega(j)^2);
            alpha1=a(j);
        end
        for k=j+1:g
            b(j)=beta1+imchi(k)*omega(k)^(2*alpha+1)/(omega(k)^2-omega(j)^2);
            beta1=b(j);
        end
    rechi(j)=2/pi*deltaomega*(a(j)+b(j))*omega(j)^(-2*alpha);
    end
    %Last element of the output: the principal part integration
    %is computed by excluding the last element of the input


    %%

    % Define the file path
filePath = 'C:\Users\Poh\Documents\ParteIm.csv';

% Read the data into a table, skipping the first 3 rows
dataTable = readtable(filePath, 'HeaderLines', 1);

%%
omega = dataTable{1:112, 1};         % eV
k      = dataTable{1:112, 2}; % segunda columna, distinto rango


chiRe=prueba(omega,k,0);

subplot(2,1,1);
x = omega;
y1 = chiRe;
plot(x,y1)

subplot(2,1,2); 
y2 = k;
plot(x,y2)
