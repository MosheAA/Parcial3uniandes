function [m_ConvMat, v_TimeArray, v_FreqTestHz] = f_GaborTFTransform(p_XIn, p_FsHz, p_F1Hz, p_F2Hz, p_FreqResHz, p_NumCycles)
%% Creado por Mario A. Valderrama Uniandes. 2020

% Transformada de Fourier método Gabor 
% Input:
%
% p_XIn  = Señal de entrada
% p_F1Hz = Límite inferior del espectrograma 
% p_F2Hz = Límite superior del espectrograma 
% p_FreqResHz = Resolución del espectrograma
% p_NumCycles = Número de ciclos de la señal de prueba

% Ejemplo:
% X = Señal. Vector (1 x tttt)
% p_F1Hz = 5
% p_F2Hz = 20
% p_FreqResHz = 0.25  
% p_NumCycles = 5


%      Creamos un vector de tiempo en segundos

v_TimeArray = 0:1/p_FsHz:length(p_XIn)/p_FsHz; 

%      Definimos un rango de frecuencias
%      las cuales usaremos para crear nuestros
%      patrones oscilatorios de prueba

v_FreqTestHz = p_F1Hz:p_FreqResHz: p_F2Hz + p_FreqResHz;

%      Creamos una matriz que usaremos para
%      almacenar el resultado de las
%      convoluciones sucesivas. En esta matriz,
%      cada fila corresponde al resultado de
%      una convolución y cada columna a todos
%      los desplazamientos de tiempo.

m_ConvMat = zeros(length(v_FreqTestHz),length(p_XIn));

%      Se obtiene la transformada de Fourier
%      de la señal p_XIn para usarla en cada iteración

p_XInfft = fft(p_XIn);

%      Ahora creamos un procedimiento iterativo
%      que recorra todas las frecuencias de prueba
%      definidas en el arreglo v_FreqTestHz

for s_FreqIter = 1:length(v_FreqTestHz)
%       Generamos una señal sinusoidal de prueba
%       que oscile a la frecuencia de la iteración
%       s_FreqIter (v_FreqTestHz[s_FreqIter]) y que tenga
%       la misma longitud que la señal p_XIn.
%       En este caso usamos una exponencial compleja.

xtest = exp(1j * 2.0 * pi * v_FreqTestHz(s_FreqIter) * v_TimeArray);

%       Creamos una ventana gaussina para
%       limitar nuestro patrón en el tiempo
%       Definimos la desviación estándar de
%       acuerdo al número de ciclos definidos
%       Dividimos entre 2 porque para un ventana
%       gaussiana, una desviación estándar
%       corresponde a la mitad del ancho de la ventana

xtestwinstd = ((1.0 / v_FreqTestHz(s_FreqIter)) * p_NumCycles) / 2.0;

%       Definimos nuestra ventana gaussiana

xtestwin = exp(-0.5 * (v_TimeArray / xtestwinstd) .^ 2.0);

%       Multiplicamos la señal patrón por
%       la ventana gaussiana
        
xtest = xtest .* xtestwin; 

%       Para cada sinusoidal de prueba obtenemos
%       el resultado de la convolución con la señal p_XIn
%       En este caso nos toca calcular la convolución
%       separadamente para la parte real e imaginaria
%       m_ConvMat[s_FreqIter, :] = np.convolve(p_XIn, np.real(xtest), 'same') + \
%                        1j * np.convolve(p_XIn, np.imag(xtest), 'same')

%       Se obtine la transformada de Fourier del patrón

fftxtest = fft(xtest);

%       Se toma únicamente la parte real para evitar
%       corrimientos de fase

fftxtest = abs(fftxtest);

%       Se obtine el resultado de la convolución realizando
%       la multiplicación de las transformadas de Fourier de
%       la señal p_XIn por la del patrón. Se obtienen las unidades de
%       entrada al cuadrado(?).
        
m_ConvMat(s_FreqIter, :) = ifft(p_XInfft .* fftxtest(1:end-1));

end
v_TimeArray = v_TimeArray - v_TimeArray(1);

end