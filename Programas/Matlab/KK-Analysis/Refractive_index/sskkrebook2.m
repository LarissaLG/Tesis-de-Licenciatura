function nreal=sskkrebook2(omega,kimag,omega1,nreal1,alpha)
    %The program inputs are 1) omega, vector of the
    % frequency (or energy) components, 2) imchi, vector of
    %the imaginary part of the susceptibility
    %under examination, 3) omega1, anchor point, 4) rechi1,
    %value of the real part at the anchor point, 5) alpha,
    %value of the moment considered.
    %The two vectors 1) and 2) must have the same length.
    %The output is the estimate of the
    %real part as obtained by using SSKK relations.
    %In order to use this program, save the whole text contained
    %in this section in a file and name it sskkrebook.m
    if size(omega,1)>size(omega,2)
        omega=omega';
    end
    if size(kimag,1)>size(kimag,2)
        nreal=nreal';
    end
    %Here the program rearranges the two vectors so that,
    %whichever their initial shape, they become row vectors.
    g=size(omega,2);
    %Size of the vectors.%
    k=0; 
    for j=1:g
        if omega(j)==omega1
            k=j;
        end
    end
    %Determination of the anchor point.
    nreal=kkrebook(omega,kimag,alpha);
    %Application of K-K relations
    nreal=nreal+omega1^(2*alpha)*omega.^(-2*alpha)*(nreal1-nreal(k));
    %The subtracted relation upgrades the estimate obtained
    %with K-K relations.