%% ASIGNATURA NEUROINGENIERIA - UNIVERSIDAD DE DE LOS ANDES
% PROGRAMA DE MAESTRÍA Y DOCTORADO EN INGENIERÍA BIOMÉDICA
% PARCIAL 3 OSCILACIONES DE ALTA FRECUENCIA (HFO)
% ESTUDIANTES:   ANDREA TORRES RUIZ - CODIGO: 200427296 & MOSHE A. AMARILLO
% PROFESOR:     PhD. MARIO ANDRES VALDERRAMA
% Artículos de referencia: Richard J. Staba, et al., J Neurophysiol 88:1743-1752, 2002
% IMPORTACION DE LA SEÑAL y SUBMUESTREO
r_Data        = load('microdata9.mat');                        % Se extraen los datos comprendiendo que estan en formato int16s_NChann     = 1;
s_FsHz        = r_Data.s_SRate;                                 % Tasa de muestreo
m_Data        = r_Data.v_MicroData;                             % Se extraen los datos en formato int16
m_Data        = double(m_Data);                                 % Convertir de int16 a double 
s_FsHz_r      = 1500;                                           % Nueva tasa de muestreo
m_Data_r      = resample(m_Data,s_FsHz_r,s_FsHz);                   % Submuestreo
v_Time        = (0:size(m_Data_r, 1) -1)./ s_FsHz_r;            % Se crea el vector de tiempo (Unidades segundos)%
clear m_Data s_FsHz r_Data
%% FILTRADO 80-500Hz
%Se diseña Filtro pasa banda con frecuencias de muestreo especificadas
cb           = 80;                                              % Límite inferior
ca           = 500;                                             % Límite superior
[b,a]        = butter(4,[cb  ca]/(s_FsHz_r/2),'bandpass');      % Se hallan coeficientes del filtro butterwoth
m_Data_F     = filtfilt(b,a,m_Data_r);                          % Se aplica el filtro a partir de los coeficientes
clear cb ca b a 
%% DETECCIÓN DE HFO (criterios: Amplitud >5SD & Duración =>6ms)
% Amplitud de ripples y Fast Ripples
window_RMS = round(s_FsHz_r*(3/1000));                          % Establecer ancho ventana deslizante de ~3ms 
m_Data_RMS = sqrt(movmean(m_Data_F.^2,window_RMS));             % Calcular RMS de la señal a través de ventana deslizante
SD_t       = 5 * std(m_Data_RMS);                               % Umbral de 5 SD 

[pks, locs] = findpeaks(m_Data_RMS,v_Time,'MinPeakDistance',...
    0.006,'MinPeakHeight', SD_t);                               % Detección de candidatos a Ripples según amplitud
% Duración de ripples 
dur_min_ms = 6;                                                 % Duración mínima de ripples 6 (milisegundos)
dur_min_sam= s_FsHz_r*(dur_min_ms/1000);                        % Duración mínima de ripples en muestras
rms_SD     = m_Data_RMS - SD_t;                                 % RMS-SD para encontrar los puntos de onset y offset del ripple
dur_max_ms = 100;
dur_max_sam= s_FsHz_r*(dur_max_ms/1000);                        % Duración máxima de ripples en muestras
min_dist_s = 10/1000;                                           % Distancia dobletes (i.e. eventos separados <10ms)
min_dist_sam = s_FsHz_r*min_dist_s;                             % Distancia dobletes (unidades muestras)
x = 1:1:length(rms_SD);
zc_RMS = ZeroX(x,rms_SD);                                       % Cruzamientos del umbral 5-SD (unidades muestras) 
zc_RMS = round(zc_RMS);                                         % Aprox cruzamientos en cero
k = 1; 
d = 1; 
for i = 1:length(locs)
    loc_temp = locs(1,i)*s_FsHz_r; % Muestras
    ini = max(zc_RMS(zc_RMS<loc_temp));
    fin = min(zc_RMS(zc_RMS>loc_temp));
        if isempty(ini) || isempty(fin)
            continue 
        end 
        if fin-ini<dur_min_sam || fin-ini>dur_max_sam           % Descartar ripples que se salgan de la duración típica
             locs_ripples_descar(1,d) = locs(i);
             pks_ripples_descar(1,d)  = pks(i);
             d = d+1;
            continue 
        end
        if k ~= 1
            dur_ripples_temp = ini-(dur_ripples(1,k-1)*s_FsHz_r); % Detectar dobletes
            if  dur_ripples_temp < min_dist_sam
                 locs_ripples_descar(1,d) = locs(i);
                 pks_ripples_descar(1,d)  = pks(i);
                 d = d+1;
                clear dur_ripples_temp
                continue 
            end
         clear dur_ripples_temp    
        end
        
    locs_ripples(1,k) = locs(i);
    pks_ripples(1,k)  = pks(i);
    time_ripples(1,k) = ini/s_FsHz_r; %tiempo de onset  (s)
    time_ripples(2,k) = fin/s_FsHz_r; %tiempo de offset (s)
    dur_ripples(1,k)  = time_ripples(2,k)-time_ripples(1,k); % Duración ripples candidatos en s
       k = k+1; 
   clear ini fin loc_temp  dur_ripples_temp
end                                       % Detección onset offset y duracción de ripples candidatos
clear k x  dur_min_ms dur_max_ms min_dist_s zc_RMS
% Ventana de observación de candidatos 
win_length_ms = 300;                                            % Ancho de ventana 300
win_length_sam= s_FsHz_r*((win_length_ms/2)/1000);              % Unidades en muestras (media ventana)
ripples_cand  = zeros(length(locs_ripples),win_length_sam*2+1 );% Preallocate
for i = 1: length(locs_ripples)
    ini = round((locs_ripples(i)*s_FsHz_r) - win_length_sam);
    fin = round((locs_ripples(i)*s_FsHz_r) + win_length_sam);
    if fin>length(m_Data_F)
        continue
    end 
    ripples_cand (i,:) = m_Data_F(ini:fin, 1);
end
%% Depuración HFOs (Contaminación 60Hz, Spikes)
% Frecuencia(s) de pico de candidatos de épocas seleccionadas
% (contaminado 60Hz)
L =  size (ripples_cand,2);
f = s_FsHz_r*(0:(L/2))/L;
freq_p = zeros(length(ripples_cand),2); % preallocate
freq_l = zeros(length(ripples_cand),2); % preallocate
for indx = 1:length(ripples_cand)
        warning off
        p_XIn       = ripples_cand (indx,:);
        if p_XIn(end) == 0
            continue
        end 
        Y           = fft(p_XIn);
        P2          = abs(Y/L);
        P1          = P2(1:L/2+1);
        P1(2:end-1) = 2*P1(2:end-1);
    [fp, fl]        = findpeaks(P1,f, 'SortStr','descend','MinPeakDistance',50);
    freq_p(indx,:)  = fp(1,1:2);
    freq_l(indx,:)  = fl(1,1:2);
    clear fp fl P1 P2 Y p_XIn
end
% figure 
% scatter(freq_l(:,1), freq_l(:,2))                             % Aquí se pueden ver eventos con picos de frecuencia en ~60Hz y sus armónicos
% Descartar ripples con contaminación de línea 60Hz y armónicos
hum_Hz = [60, 120, 180, 240, 300, 360, 420]; 
g = 1; 
for i = 1:length(locs_ripples)
         if freq_l(i,1) == 0
           continue 
         end 
          for l=1:length(hum_Hz)
            ini = hum_Hz(l)-3;
            fin = hum_Hz(l)+3;
            cond(l) = (freq_l(i,1)>ini && freq_l(i,1)<fin)|| (freq_l(i,2)>ini && freq_l(i,2)<fin);
          end 
          if all(~cond) % Si el primer o segundo pico de frecuencia  no está entre el ruido de línea y sus armónicos...                 
                    locs_ripples_corr(1,g) = locs_ripples(1,i);             % Localización del pico (s)
                    pks_ripples_corr(1,g)  = pks_ripples(1,i);              % Amplitud del pico ()
                    time_ripples_corr(1,g) = time_ripples(1,i);             % Time onset  (s)
                    time_ripples_corr(2,g) = time_ripples(2,i);             % Time offset (s)
                    dur_ripples_corr(1,g)  = dur_ripples(1,i) ;      
                    freq_p_1(g,:)  = freq_p(i,1:2);
                    freq_l_1(g,:)  = freq_l(i,1:2); 
                    g = g+1;
          else
             locs_ripples_descar(1,d) = locs_ripples(i);
             pks_ripples_descar(1,d)  = pks_ripples(i);
             d = d+1;
              continue 
          end 
  clear cond ini fin
end
clear indx freq_p freq_l cond l g hum_Hz p_XIn freq_p freq_l
% Épocas depuradas
clear ripples_cand
ripples_cand  = zeros(length(locs_ripples_corr),win_length_sam*2+1 );% Preallocate
for i = 1: length(locs_ripples_corr)
    ini = round((locs_ripples_corr(i)*s_FsHz_r) - win_length_sam);
    fin = round((locs_ripples_corr(i)*s_FsHz_r) + win_length_sam);
    if fin>length(m_Data_F)
        continue
    end 
    ripples_cand (i,:) = m_Data_F(ini:fin, 1);
    clear ini fin
end
wideband_cand = zeros(length(locs_ripples_corr),win_length_sam*2+1 );% Preallocate
for i = 1: length(locs_ripples_corr)
ini = round((locs_ripples_corr(i)*s_FsHz_r) - win_length_sam);
fin = round((locs_ripples_corr(i)*s_FsHz_r) + win_length_sam);
if fin>length(m_Data_F)
continue
end
wideband_cand (i,:) = m_Data_r(ini:fin, 1);
    clear ini fin
end
%% Descartar ripples que tengan ruido por spike en wideband
diff_wideband = zeros(size(locs_ripples_corr,2),1);             % preallocate
for i = 1: size(locs_ripples_corr,2)
    diff_wideband(i,:) = max(diff(wideband_cand(i,:)));
end                         % Valor máximo de la derivada wideband
freq_p_wideband = zeros(size(locs_ripples_corr,2),2);           % preallocate
freq_l_wideband = zeros(size(locs_ripples_corr,2),2);           % preallocate   
for indx = 1:size(locs_ripples_corr,2)
        warning off
        p_XIn       = wideband_cand (indx,:);
        if p_XIn(end) == 0
            continue
        end 
        Y           = fft(p_XIn);
        P2          = abs(Y/L);
        P1          = P2(1:L/2+1);
        P1(2:end-1) = 2*P1(2:end-1);
    [fp, fl]        = findpeaks(P1,f, 'SortStr','descend');
    freq_p_wideband(indx,:)  = fp(1,1:2);
    freq_l_wideband(indx,:)  = fl(1,1:2);
    clear fp fl P1 P2 Y p_XIn
end                       % Picos de frecuencia dominante Wideband 

Feat = cat (2,diff_wideband,freq_p_wideband(:,1),...
    freq_p_wideband(:,2),freq_l_wideband(:,1),...
    freq_l_wideband(:,2));                                      % Características para incluir en PCA

[~, score]= pca(Feat);                                          % PCA
comp = [1, 2];                                                  % Componentes elegidos 
cat_comp = cat(2, score(:,comp(1)), score(:,comp(2)));          % Concatenar componentes elegidos
Z = linkage(cat_comp,'centroid', 'cityblock'	);
clus = cluster(Z,'Maxclust',8);
% figure
% scatter(cat_comp(:,1),cat_comp(:,2),10,clus)
k = mode (clus);     %cluster con más eventos (en teoría el que más eventos válidos tiene
for i = 1:max(clus) 
    clusters{i,:} = find(clus==i);                               
end  % índices de los eventos
for i = k:k         
    clus_temp = clusters{i,1};
scatter(score(clus_temp,comp(1)), ...
    score(clus_temp,comp(2)))
hold on 
clear clus_temp
end  % Plot cluster 1
for i = 1:max(clus) 
    if i== k
        continue
    end 
        clus_temp = clusters{i,1};
    scatter(score(clus_temp,comp(1)), ...
        score(clus_temp,comp(2)), 'x', 'k')
    hold on 
    clear clus_temp
        xlabel(['PC',num2str(comp(1))]),...
        ylabel(['PC',num2str(comp(2))]),...
        grid on

clear clus_temp
end  % Plot clusters restantes

% for i = 1:max(clus)
%    clus_temp = clusters{i,1};
%    figure
%    for c = 1:round(length(clus_temp)/2)
%    subplot(121)
%    plot(wideband_cand (clus_temp(c,1),:))
%    grid on
%    subplot(122)
%    plot(ripples_cand (clus_temp(c,1),:))
%    title(['Ejemplos de HFOs incluidos en el cluster:', num2str(i)])
%    grid on
%    hold on
%  pause
%    end 
%    pause
% clear clus_temp
% end  % Plot de visualización de algunos eventos por cluster
% Actualizando 
for l = 1:clus 
    clus_temp = clusters{l,1};
    if l == k
            locs_ripples_fin(1,:) = locs_ripples_corr(1,clus_temp);             % Localización del pico (s)
            pks_ripples_fin(1,:)  = pks_ripples_corr(1,clus_temp);              % Amplitud del pico ()
            time_ripples_fin(1,:) = time_ripples_corr(1,clus_temp);             % Time onset  (s)
            time_ripples_fin(2,:) = time_ripples_corr(2,clus_temp);             % Time offset (s)
            dur_ripples_fin(1,:)  = dur_ripples_corr(1,clus_temp) ;      
    else
        for s = 1:size(clus_temp,1)
            locs_ripples_descar(1,d) = locs_ripples_corr(clus_temp(s));
            pks_ripples_descar(1,d)  = pks_ripples_corr(clus_temp(s));
            d = d+1;
        end 
    end 
    clear clus_temp
end 
clear ripples_cand wideband_cand
ripples_cand  = zeros(size(locs_ripples_fin,2),win_length_sam*2+1);% Preallocate
for i = 1: size(locs_ripples_fin,2)
    ini = round((locs_ripples_fin(i)*s_FsHz_r) - win_length_sam);
    fin = round((locs_ripples_fin(i)*s_FsHz_r) + win_length_sam);
    if fin>length(m_Data_F)
        continue
    end 
    ripples_cand (i,:) = m_Data_F(ini:fin, 1);
    clear ini fin
end
wideband_cand = zeros(size(locs_ripples_fin,2),win_length_sam*2+1);% Preallocate
for i = 1: size(locs_ripples_fin,2)
    ini = round((locs_ripples_fin(i)*s_FsHz_r) - win_length_sam);
    fin = round((locs_ripples_fin(i)*s_FsHz_r) + win_length_sam);
    if fin>length(m_Data_F)
        continue
    end
wideband_cand (i,:) = m_Data_r(ini:fin, 1);
    clear ini fin
end
RMS_cand      = zeros(size(locs_ripples_fin,2),win_length_sam*2+1);% Preallocate
for i = 1: size(locs_ripples_fin,2)
    ini = round((locs_ripples_fin(i)*s_FsHz_r) - win_length_sam);
    fin = round((locs_ripples_fin(i)*s_FsHz_r) + win_length_sam);
    if fin>length(m_Data_F)
        continue
    end
RMS_cand (i,:) = m_Data_RMS(ini:fin, 1);
    clear ini fin
end
    L =  size (ripples_cand,2);                                 % Calcular la fft de los eventos elegidos
    f = s_FsHz_r*(0:(L/2))/L;
    freq_p = zeros(length(ripples_cand),2); % preallocate
    freq_l = zeros(length(ripples_cand),2); % preallocate
for i = 1: size(locs_ripples_fin,2)
        warning off
        p_XIn       = ripples_cand (i,:);
        if p_XIn(end) == 0
            continue
        end 
        Y           = fft(p_XIn);
        P2          = abs(Y/L);
        P1          = P2(1:L/2+1);
        P1(2:end-1) = 2*P1(2:end-1);
    [fp, fl]        = findpeaks(P1,f, 'SortStr','descend','MinPeakDistance',50);
    freq_p(i,:)  = fp(1,1:2);
    freq_l(i,:)  = fl(1,1:2);
    clear fp fl P1 P2 Y p_XIn
end


%% División Ripples y FR

R  = find(freq_l(:,1)>80 & freq_l(:,1)<180 | freq_l(:,2)>80);
FR = find(freq_l(:,1)>180 & freq_l(:,1)<400 & freq_l(:,2)>150);

Dur_Ripl     = dur_ripples_fin(1, R);
Dur_F_Ripl   = dur_ripples_fin(1, FR);

freq_Ripl     = freq_l(R,1);
freq_F_Ripl   = freq_l(FR,1);
freq_peak_Ripl     = freq_p(R,1);
freq_peak_F_Ripl   = freq_p(FR,1);

Peak_all = cat(1,freq_peak_Ripl,freq_peak_F_Ripl);
Freq_all = cat(1,freq_Ripl,freq_F_Ripl);
Dur_all  =  cat(2,Dur_Ripl,Dur_F_Ripl);
%% Verificación visual. 1 por 1 de los HFOs que fueron detectados 
for indx = 70:170
 %indx = 1;% length(ripples_cand)
figure 
subplot(311)
plot(RMS_cand (indx,:))              % Plot RMS
grid on

subplot(312)
plot(wideband_cand (indx,:))         % Plot Wideband
grid on

subplot(313)
plot(ripples_cand (indx,:))           % Plot Ripples
grid on
pause
close all

end 
%% Plots de ripples válidos 
% Parámetros espectrograma de ripples
for i =8
indx        = R(i);
p_XIn       = ripples_cand (indx,:);
p_FsHz      = s_FsHz_r;
p_F1Hz      = 60;
p_F2Hz      = 500;
p_FreqResHz = 1;  
p_NumCycles = 5;
% Espectrograma método Gabor
[m_ConvMat, ~, v_FreqTestHz] = ...
    f_GaborTFTransform(p_XIn, p_FsHz, p_F1Hz, p_F2Hz, ...
    p_FreqResHz, p_NumCycles);                                  % Espectrograma 

v_TimeArray = (-win_length_ms/2):(1/p_FsHz)*1000:(win_length_ms/2);
% Valores proyección del espectrograma en dominio de la frecuencia 
 for i = 1:size(m_ConvMat,1)
    p_fft(i) =max(abs(m_ConvMat(i,:)));
 end 
figure 
surf(v_TimeArray,v_FreqTestHz,abs(m_ConvMat),...       % Plot de superficie
    'EdgeColor','none');   
 axis xy; axis tight; colormap(jet); view(0,90);               % Propiedades de visualización
    xlabel('Tiempo (ms)'),...
    ylabel('Frecuencia (Hz)'),...
    zlabel('Energía (U.A.)')
    zlim([0 1500])
 % Proyección dominio tiempo 
 hold on
plot3(v_TimeArray(1:size(ripples_cand (indx,:),2)),500*ones(1,size(ripples_cand (indx,:),2)), (ripples_cand (indx,:)*2.5)+500, 'LineWidth',1, 'Color' , 'k')
% Proyección dominio frecuencia
 hold on
plot3(min(v_TimeArray)*ones(1,size(v_FreqTestHz,2)),v_FreqTestHz,p_fft, 'LineWidth',1, 'Color' , 'k')
%pause 
%close all
end 
%exportgraphics(gcf,'TF-gaborBIaplastado.emf','Resolution',1000)
%%% En paper: Bimodal:R(8) 
%             Unimodal baja frecuencia: R(9)
%             Unimodal alta frecuencia: FR(5)
% Para ver las figuras, cambie el valor de i en la línea 347 y R o FR en la línea 348 según corresponda 
