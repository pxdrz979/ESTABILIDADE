clc, clear, close;
%% Puxar dados de planilhas
[file1, path1] = uigetfile('*.txt','Selecione Arquivos com dados do XFLR5');
if isequal(file1, 0)
    error('Arquivo de dados não selecionado');
end
arquivo1 = fullfile(path1, file1);

[file2, path2] = uigetfile('*.txt','Selecione arquivo com os dados CFD');
if isequal(file2, 0)
    error('Arquivo de dados não selecionado');
end
arquivo2 = fullfile(path2, file2);
    [alpha1cpt, CLcpt, CDcpt, Cmcpt, alpha1w, Lw, Dw, Mw] = importaDadosAsa(arquivo1, arquivo2);

[file3, path3] = uigetfile('.txt', 'Selecione arquivo com os dados do EH');
if isequal(file3, 0)
    error('Arquivo de dados não selecionado');
end 
arquivo3 = fullfile(path3, file3);
    [alpha1eh, Leh, Deh, Meh] = importaDadosEh(arquivo3);

%% Dados Fisicos
% Dados do ar
dens_ar = 1.180724; % densidade do ar [kg/m³]
visc_ar = 0.000018325; % viscosidade do ar [Pa.s]
vel = 15; % velocidade de voo [m/s]
kvisc = visc_ar/dens_ar; % viscosidade cinematica do ar
% Dados da asa
cw = 0.2909063; % Corda média aerodinamica [m]
Sw = 0.934625; % Área da Asa [m²]
bw = 2.25; % Envergadura da asa [m]
ar_asa = (bw^2)/Sw; % Razão de alongamento da asa
e_obert_asa = 1/(1.05+0.007*pi*ar_asa); % eficiencia da asa
kcdi_asa = 1/(pi*e_obert_asa*ar_asa);
ac_asa = 0.25*cw;
cp_asa = 2/3*cw;
% Dados Empenagem Horizontal
Seh = 0.22756; % Área da EH [m²]
beh = 0.64; % Envergadura da EH [m]
ceh = Seh/beh; % Corda média aerodinamica EH [m]
ar_eh = (beh^2)/Seh; % Razão de alongamento do EH
e_obert_eh = 1/(1.05+0.007*pi*ar_eh); % eficiencia do EH
kcdi_eh = 1/(pi*e_obert_eh*ar_eh);
ac_eh = 0.25*ceh; % Calculo do centro aerodinamico
cp_eh = 2/3*ceh;
%% Configurações da Simulação
% Tamanho de Dominio 
H_Asa = cw*1000; % Comprimento Corda para tamanho de dominio [mm]
W_Asa = bw*1000; % Envergadura de Asa [mm]
L_Domain = 5*H_Asa + 10*H_Asa; % Comprimento do dominio [mm]
H_Domain = 5*H_Asa; % Altura do dominio [mm]
W_Domain = 10*H_Asa; % Largura do dominio [mm]
% Velocidade
deginrad = deg2rad(5);
vely = vel*sin(deginrad);
velx = vel*cos(deginrad);
cordy = sin(deginrad);
cordx = cos(deginrad);
% Calculo para Camada Limite
rey = (dens_ar*vel*ceh)/visc_ar; % Reynolds do escoamento
yplus = 1;
cf = (2*log10(rey)-0.65)^(-2.3); % Coeficiente de atrito da parede
tw = 0.5*cf*dens_ar*(vel^2); % Atrito da parede
ut = sqrt(tw/dens_ar); % Velocidade de fricção
y1 = ((yplus*visc_ar)/(ut*dens_ar)); % Altura da primeira camada
delta = 0.14*(kvisc/vel)*(rey/log(rey))*1.5*log(rey); % Boundary Layer Thickness
n_cam = lognabase12(1-((delta/y1)*(1-1.2)))-1; % Numero de Camadas do Inflation
y1mm = y1*1000;
%% DADOS X5
% Converte dados em table e define parametros
table_x5 = [array2table(alpha1cpt) array2table(CLcpt) array2table(CDcpt) array2table(Cmcpt)];
alpha0x5 = table_x5.alpha1cpt==0; % Define linha onde o angulo é 0
alphastolrowx5 = table_x5.alpha1cpt==14; % Define linha onde o angulo é onde há CL maximo
% Dados importantes a serem obtidos
CLx5max = table_x5.CLcpt(alphastolrowx5);
CDx5max = table_x5.CDcpt(alphastolrowx5);
CL0x5 = table_x5.CLcpt(alpha0x5);
CD0x5 = table_x5.CDcpt(alpha0x5);
%% Dados CFD ASA
% Calculo dos Coeficientes
CDw = Dw./(0.5*(vel^2)*Sw*dens_ar);
CLw = Lw/(0.5*(vel^2)*Sw*dens_ar);
CMw = Mw/(0.5*(vel^2)*Sw*dens_ar*cw);
mac_w = -14.6783348;
cmac_w = mac_w/(0.5*(vel^2)*Sw*dens_ar*cw);
% Converte dados em table e define parametros
table_cfd_w = [array2table(alpha1w) array2table(CLw) array2table(CDw) array2table(CMw) array2table(Lw) array2table(Dw) array2table(Mw)];
alpha0w = table_cfd_w.alpha1w==0; % Define linha onde o angulo é 0
alphastolrow_w = table_cfd_w.CLw==max(table_cfd_w.CLw); % Define linha onde o angulo é onde há CL maximo
% Loop substituindo NaN por 0
for varName = table_cfd_w.Properties.VariableNames
    % Get the column data
    columnData = table_cfd_w.(varName{1});
    
    % Check if the column is numeric
    if isnumeric(columnData)
        % Replace NaN with 0
        columnData(isnan(columnData)) = 0;
        % Assign the modified column back to the table
        table_cfd_w.(varName{1}) = columnData;
    end
    % For non-numeric data, you can define other replacements if necessary
end
% Dados importantes a serem obtidos
alphaclw = (table_cfd_w.CLw(table_cfd_w.alpha1w==5)-table_cfd_w.CLw(table_cfd_w.alpha1w==1))/(5-1); 
alphacmw = (table_cfd_w.CMw(table_cfd_w.alpha1w==5)-table_cfd_w.CMw(table_cfd_w.alpha1w==1))/(5-1);
alphastol_w = table_cfd_w.alpha1w(alphastolrow_w); % Angulo de estol da asa
CLmaxw = table_cfd_w.CLw(alphastolrow_w); % CL maximo
CDmaxw = table_cfd_w.CDw(alphastolrow_w); % CD maximo
CM0w = table_cfd_w.CMw(alpha0w); % CM a angulo 0
CL0w = table_cfd_w.CLw(alpha0w); % CL a angulo 0
CD0w = table_cfd_w.CDw(alpha0w); % CD a angulo 0
%% Dados CFD EH
% Calculo dos Coeficientes
CDeh = Deh./(0.5*(vel^2)*Seh*dens_ar);
CLeh = Leh/(0.5*(vel^2)*Seh*dens_ar);
CMeh = Meh/(0.5*(vel^2)*Seh*dens_ar*ceh);
mac_eh = -14.6783348;
cmac_eh = mac_eh/(0.5*(vel^2)*Seh*dens_ar*cw);
% Converte dados em table e define parametros
table_cfd_eh = [array2table(alpha1eh) array2table(CLeh) array2table(CDeh) array2table(CMeh) array2table(Leh) array2table(Deh) array2table(Meh)];
alpha0eh = table_cfd_eh.alpha1eh==0; % Define linha onde o angulo é 0
alphastolrow_eh = table_cfd_eh.CLeh==max(table_cfd_eh.CLeh); % Define linha onde o angulo é onde há CL maximo
% Loop substituindo NaN por 0
for varName = table_cfd_eh.Properties.VariableNames
    % Get the column data
    columnData = table_cfd_eh.(varName{1});
    
    % Check if the column is numeric
    if isnumeric(columnData)
        % Replace NaN with 0
        columnData(isnan(columnData)) = 0;
        % Assign the modified column back to the table
        table_cfd_eh.(varName{1}) = columnData;
    end
    % For non-numeric data, you can define other replacements if necessary
end
% Dados importantes a serem obtidos
alphacleh = (table_cfd_eh.CLeh(table_cfd_eh.alpha1eh==5)-table_cfd_eh.CLeh(alpha0eh))/(5-0); 
alphacmeh = (table_cfd_eh.CMeh(table_cfd_eh.alpha1eh==5)-table_cfd_eh.CMeh(alpha0eh))/(5-0);
alphastol_eh = table_cfd_eh.alpha1eh(alphastolrow_eh); % Angulo de estol da asa
CLmaxeh = table_cfd_eh.CLeh(alphastolrow_eh); % CL maximo
CDmaxeh = table_cfd_eh.CDeh(alphastolrow_eh); % CD maximo
CM0eh = table_cfd_eh.CMeh(alpha0eh); % CM a angulo 0
CL0eh = table_cfd_eh.CLeh(alpha0eh); % CL a angulo 0
CD0eh = table_cfd_eh.CDeh(alpha0eh); % CD a angulo 0
%% Print Resultados
clc;
fprintf('=======================================\nDados Fisico Asa \n')
fprintf('\nCoeficiente de Proporcionalidade = %.5f\n',kcdi_asa)
fprintf('Coeficiente de Oswald = %.5f\n',e_obert_asa)
fprintf('Razão de Alongamento = %.5f\n',ar_asa)
fprintf('\n=======================================\nDominio do seu enclosure \n')
fprintf('Comprimento do dominio = %.7f mm\n',L_Domain)
fprintf('Altura do dominio = %.7f mm\n',H_Domain)
fprintf('Largura do dominio = %.7f mm\n',W_Domain)
fprintf('\nDados para Inflation\n')
fprintf('Altura da primeira camada = %.7f mm\n',y1mm)
fprintf('Numero de Camadas do Inflation = %.0f\n',n_cam)
fprintf('Velocidade em Y = %.7f\n',vely)
fprintf('Velocidade em X = %.7f\n',velx)
fprintf('Valor em Y = %.7f\n',cordy)
fprintf('Valor em X = %.7f\n',cordx)
fprintf('\n=======================================\nDados da Asa X5\n')
fprintf('\nCoeficiente de Sustentação a angulo 0(CL0) = %.7f\n',CL0x5)
fprintf('Coeficiente de Arrasto a angulo 0(CD0) = %.7f\n',CD0x5)
fprintf('Coeficiente de Sustentação Maximo(CLmax) = %.7f\n',CLx5max)
fprintf('Coeficiente de Arrasto Maximo(CDmax) = %.7f\n',CDx5max)
fprintf('\n=======================================\nDados da Asa CFD\n')
fprintf('\nAngulo de Estol da Asa = %.0f\n',alphastol_w)
fprintf('Coeficiente Angular da Asa = %.7f\n',alphaclw)
fprintf('Coeficiente Angular de Momento da Asa = %.7f\n',alphacmw)
fprintf('Coeficiente de Sustentação a angulo 0(CL0) = %.7f\n',CL0w)
fprintf('Coeficiente de Arrasto a angulo 0(CD0) = %.7f\n',CD0w)
fprintf('Coeficiente de Sustentação Maximo(CLmax) = %.7f\n',CLmaxw)
fprintf('Coeficiente de Arrasto Maximo(CDmax) = %.7f\n',CDmaxw)
fprintf('Coeficiente de Momento(CM0) = %.7f\n',CM0w)
fprintf('Coeficiente de Momento no Centro Aerodinamico = %.7f\n',cmac_w)
fprintf('\n=======================================\nDados do EH CFD\n')
fprintf('\nAngulo de Estol do EH = %.0f\n',alphastol_eh)
fprintf('Coeficiente Angular do EH = %.7f\n',alphacleh)
fprintf('Coeficiente Angular de Momento do EH = %.7f\n',alphacmeh)
fprintf('Coeficiente de Sustentação a angulo 0(CL0) = %.7f\n',CL0eh)
fprintf('Coeficiente de Arrasto a angulo 0(CD0) = %.7f\n',CD0eh)
fprintf('Coeficiente de Sustentação Maximo(CLmax) = %.7f\n',CLmaxeh)
fprintf('Coeficiente de Arrasto Maximo(CDmax) = %.7f\n',CDmaxeh)
fprintf('Coeficiente de Momento(CM0) = %.7f\n',CM0eh)
fprintf('Coeficiente de Momento no Centro Aerodinamico = %.7f\n',cmac_eh)

%% Graficos
close;
tiledlayout(2,2)
nexttile
plot(alpha1w,CLw)
title('CL x alpha')
xlabel('alpha')
ylabel('CL')
nexttile
plot(alpha1w,CDw)
title('CD x alpha')
xlabel('alpha')
ylabel('CD')
nexttile
plot(alpha1w,CMw)
title('CM x alpha')
xlabel('alpha')
ylabel('CM')
nexttile
plot(CDw, CLw)
title('Cl x Cd')
xlabel('CD')
ylabel('CL')
%% Funções
function [alpha1cpt, CLcpt, CDcpt, Cmcpt, alpha1w, Lw, Dw, Mw] = importaDadosAsa(arquivo1, arquivo2)
    %Importa dados da Asa do XFLR5
    dadoscpt = readtable(arquivo1); %Importa dados do txt como uma tabela
    dadoscfd = readtable(arquivo2); %Importa dados do xlsx como uma tabela
    %Extrai colunas de interesse e transforma em arrays numericos
    alpha1cpt = table2array(dadoscpt(:,"alpha"));
    CLcpt = table2array(dadoscpt(:,"CL"));
    CDcpt = table2array(dadoscpt(:,"CD"));
    Cmcpt = table2array(dadoscpt(:,"Cm"));
    alpha1w = table2array(dadoscfd(:,"alpha"));
    Lw = table2array(dadoscfd(:,"L"));
    Dw = table2array(dadoscfd(:,"D"));
    Mw = table2array(dadoscfd(:,"M"));
end
function [alpha1eh, Leh, Deh, Meh] = importaDadosEh(arquivo3)
    dadoseh = readtable(arquivo3);
    alpha1eh = table2array(dadoseh(:,"alpha"));
    Leh = table2array(dadoseh(:,"L"));
    Deh = table2array(dadoseh(:,"D"));
    Meh = table2array(dadoseh(:,"M"));
end
function [log12] = lognabase12(x)

log12 = log(x)./log(1.2);

end

