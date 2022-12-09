suitable_init_val=@(x) blasius(x)-1;
temp=fzero(suitable_init_val,[0,1]); % to find value of f2(0) for which f1(∞) is 1

fprintf('Optimized Solution is %f\n',temp);
    
function result=blasius(x)
eta=linspace(0,7,1000000);
n=length(eta);
f=ones(n,3); % matrix of f0, f1, and f2
f(1,1)=0;f(1,2)=0;f(1,3)=x; 
der={@(f0,f1,f2) f1; @(f0,f1,f2) f2; @(f0,f1,f2) -0.5*f0*f2}; % cell of functions 
h=eta(2)-eta(1); % step size
%% Runge Kutta method
for i=1:n-1
    for j=1:3
        d=der{j};
        k1=d(f(i,1),f(i,2),f(i,3));
        k2=d(f(i,1)+1/2*k1*h,f(i,2)+1/2*k1*h,f(i,3)+1/2*k1*h);
        k3=d(f(i,1)+1/2*k2*h,f(i,2)+1/2*k2*h,f(i,3)+1/2*k2*h);
        k4=d(f(i,1)+k3*h,f(i,2)+k3*h,f(i,3)+k3*h);
        f(i+1,j)=f(i,j)+1/6*h*(k1+2*k2+2*k3+k4);
    end
end
result=f(end,2);
%% Plotting results
figure(1);
grid('on'); % to plot f0, f1, and f2 in the same grid
xlim([0 2]);
plot(f(:,1),eta,'r-','LineWidth',2);
hold('on');plot(f(:,2),eta,'b-','LineWidth',2);
hold('on');plot(f(:,3),eta,'k-','LineWidth',2);
xlabel('f, f'' and f''''', 'FontSize', 15);
ylabel('\eta', 'FontSize', 15);
legend('f(\eta)', 'f''(\eta)', 'f''''(\eta)', 'Location','best');
hold('off'); % to remove the previous f0, f1, f2 curves
legend boxoff 
set(gca,'FontName','Times New Roman','Fontsize',18,'XColor','k','YColor','k','LineWidth',1);
end    