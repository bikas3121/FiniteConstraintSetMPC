% plot simulation results from python 
% read data
data = readmatrix('UF-4_100_1000_10000_23Nov30_14-37.csv'); 


t = data(:,1);  % simulation time samples


%% Unfiltered Data
% Reference signal
ref = data(:,2);   % filtered reference signal

% Direct quantization
u_direct = data(:,3);  % directly quantized reference signal with unifrom quantizer
u_directINL = data(:,4);  % directly quantized reference signal with non-unifrom quantizer

% Closedform ideal uniform MPC
u_mhoq = data(:,5);  % directly quantized reference signal with non-unifrom quantizer
u_mhoq_INL = data(:,6);  % directly quantized reference signal with non-unifrom quantizer

% Ideal uniform MPC 
u_mpc = data(:,7); % optimally quantized reference signal with unifrom quantizer without INL feedback
u_mpc_INL = data(:,8); % optimally quantized reference signal with non-unifrom quantizer without INL feedback

% Nonideal uniform MPC
nu_mpcINL = data(:,9); % optimally quantized reference signal with non-unifrom quantizer with INL feedback
nu_mpcINL_trun = data(:,10); % optimally quantized reference signal with non-unifrom quantizer with INL feedback

%%  Filter and Filtered signal 
Fc = 1e2;
Fs = 1e3;
Wn = Fc/(Fs/2);
[b,a] = butter(2,Wn, 'low');

filt_ref = filter(b,a, ref);

% Direct quantization
filt_u_direct = filter(b,a,u_direct);
filt_u_directINL = filter(b,a,u_directINL);

filt_u_mhoq = filter(b,a,u_mhoq);
filt_u_mhoq_INL = filter(b,a,u_mhoq_INL);

filt_u_mpc = filter(b,a,u_mpc);
filt_u_mpc_INL = filter(b,a,u_mpc_INL);

filt_nu_mpcINL = filter(b,a,nu_mpcINL);
filt_nu_mpcINL_trun = filter(b,a,nu_mpcINL_trun);

%% Error and variance
err_direct = filt_ref-filt_u_direct;
var_direct = var(err_direct);

err_directINL = filt_ref-filt_u_directINL;
var_directINL = var(err_directINL);

err_mhoq = filt_ref - filt_u_mhoq;
var_mhoq = var(err_mhoq);

err_mhoqINL = filt_ref - filt_u_mhoq_INL;
var_mhoqINL = var(err_mhoqINL);
 
err_mpc = filt_ref - filt_u_mpc;
var_mpc = var(err_mpc);

err_mpcINL = filt_ref - filt_u_mpc_INL;
var_mpcINL = var(err_mpcINL);
 
err_nu_mpcINL = filt_ref - filt_nu_mpcINL;
var_nu_mpcINL = var(err_nu_mpcINL);

err_nu_mpcTrun = filt_ref - filt_nu_mpcINL_trun;
var_nu_mpcINLTrun = var(err_nu_mpcTrun);

% x = categorical({'D','D-INL',  'cf-MPC','cf-MPC-INL','U-MPC','U-MPC-INL','NU-MPC-INL','NU-MPC-TrunINL'});
x = categorical({'Direct','Cf-MPC','MPC','MPC-INL'});
x = reordercats(x,{'Direct','Cf-MPC','MPC','MPC-INL'});
y = [round(var_direct,4) round(var_directINL,4); round(var_mhoq,4) round(var_mhoqINL,4);  round(var_mpc,4) round(var_mpcINL,4); round(var_nu_mpcINL,4) round(var_nu_mpcINLTrun,4) ];
figure()
% b = bar(x,y);
b = bar(x,y,'FaceColor','flat');
b(1).CData(1,:) = [0 0.4470 0.7410];
b(1).CData(2,:) = [0 0.4470 0.7410];
b(1).CData(3,:) = [0 0.4470 0.7410];
b(1).CData(4,:) = [0.8500 0.3250 0.0980];
b(2).CData = [0.8500 0.3250 0.0980];
xtips1 = b(1).XEndPoints;
ytips1 = b(1).YEndPoints;
labels1 = string(b(1).YData);
text(xtips1,ytips1,labels1,'HorizontalAlignment','center','VerticalAlignment','bottom')
xtips2 = b(2).XEndPoints;
ytips2 = b(2).YEndPoints;
labels2 = string(b(2).YData);
text(xtips2,ytips2,labels2,'HorizontalAlignment','center','VerticalAlignment','bottom')
xtickangle(45)
ylabel('Error Variance')
legend('Ideal DAC (Uniform)', 'Non-ideal DAC (INL)')
grid minor

%%
% x = categorical({'D','D-INL',  'cf-MPC','cf-MPC-INL','U-MPC','U-MPC-INL','NU-MPC-INL','NU-MPC-TrunINL'});
% x = categorical({'Direct','Cf-MPC','MPC','MPC-INL'});
% x = reordercats(x,{'Direct','Cf-MPC','MPC','MPC-INL'});
% y = [round(var_direct,4) round(var_directINL,4); round(var_mhoq,4) round(var_mhoqINL,4);  round(var_mpc,4) round(var_mpcINL,4); round(var_nu_mpcINL,4) round(var_nu_mpcINLTrun,4) ];
% figure()
% b = bar(y,'FaceColor','flat');
% b(1).CData(1,:) = [0 0.4470 0.7410];
% b(1).CData(2,:) = [0 0.4470 0.7410];
% b(1).CData(3,:) = [0 0.4470 0.7410];
% b(1).CData(4,:) = [0.8500 0.3250 0.0980];
% b(2).CData = [0.8500 0.3250 0.0980];
% 

