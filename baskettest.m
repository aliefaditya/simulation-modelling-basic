clear all;
close all;
clc;

%inisialisasi nilai awal

g = -9.8; % gravitasi (m/s^2)
x(1) = 0; % titik awal x (numerik)
y(1) = 1.7; % titik awal y = tinggi pemain (numerik)
x2(1) = 0; % titik awal x2 (analitik)
y2(1) = 1.7; % titik awal y2 = tinggi pemain (analitik)
degree = 37; % derajat
v0 = 12.5; % kecepatan awal (m/s)
vx = v0*cos(degree* pi/180); % kecepatan x == deraajat dirubah ke radian
vy(1) = v0*sin(degree* pi/180); % kecepatan y
tAkhir = (-2*vy(1))/g; %waktu sampai parabola sempurna 
Dt = 0.01;
T = 0:Dt:tAkhir;

for i=2:columns(T)
  
  vy(i) = vy(i-1) + (g*Dt); %kecepatan sebelumnya + dipengaruhi gaya gravitasi
  x(i) = x(i-1) + (vx*Dt); %posisi tergantung kecepatan terhadap waktu, vx selalu konstan
  y(i) = y(i-1) + (vy(i) * Dt); 
  if((x(i) >= 11.5) && (x(i) <= 12.5) && (y(i) >= 3.5) && (y(i) <= 3.7))
    disp("Bola menabrak ring basket (solusi Numerik)");
    yjatuh(1) = y(i);
    vjatuh(1) = 0;
    tAkhirjatuh = 0.86851; %(2h/g)^0.5
    Tjatuh = 0:Dt:tAkhirjatuh;
    vjatuh(2) = vjatuh(1) + g*Dt;
    yjatuh(2) = yjatuh(1) + vjatuh(2)*Dt;
    for j=3:columns(Tjatuh)
      yjatuh(j) = g*(Dt^2) + 2*yjatuh(j-1) - yjatuh(j-2);
    end
    break
  end
end

%solusi analitik

for m=2:columns(T)
  x2(m) = vx*T(m); 
  y2(m) = 1/2*g*(T(m)^2) + v0*sin(degree*pi/180)*T(m) + y2(1);
  if((x2(m) >= 11.5) && (x2(m) <= 12.5) && (y2(m) >= 3.5) && (y2(m) <= 3.7))
    disp("Bola menabrak ring basket (solusi Analitik)");
    yjatuh2(1) = y2(m); %ketinggian
    vjatuh2 = 0;
    tAkhir2 = 0.86558; %waktu jatuh
    Tjatuh2 = 0:Dt:tAkhir2;
    for n=2:columns(Tjatuh2)
      yjatuh2(n) = 0 + 1/2*g*(Tjatuh2(n)^2) + yjatuh2(1);
    end
    break
  end
end

figure();
plot(0,0);
plot(x,y); % Gambar lintasan
hold on;
plot(11.7,3.7,'b','markersize',70);
p = plot(x(1),y(1),'r.','markersize', 30); % buat objek bola
q = plot(x2(1),y2(1),'g.','markersize', 30); % buat objek bola
plot(x(i),y2(m));
pause(Dt);

for k=2 : length(y2)
   if (k<=119)
     delete(p);
     p = plot(x(k),y(k),'r.','markersize', 30); % buat objek bola
   end
   delete(q);
   q = plot(x2(k),y2(k),'g.','markersize', 30); % buat objek bola
   pause(Dt);
end

delete(p);
delete(q);

for l=1 : length(yjatuh)
  r = plot(x(i),yjatuh(l),'r.','markersize',30);
  s = plot(x2(m),yjatuh2(l),'g.','markersize', 30); % buat objek bola
  pause(Dt);
  delete(r);
  delete(s);
end


r = plot(x(i),yjatuh(j),'r.','markersize',30);
s = plot(x2(m),yjatuh2(n),'g.','markersize', 30); % buat objek bola

disp('Jarak tempuh x solusi Numerik =');
disp(x(i));
disp('Jarak tempuh y solusi Numerik =');
disp(y(i)); 
disp('---------------------------------');
disp('Jarak tempuh x solusi Analitik =');
disp(x2(m));
disp('Jarak tempuh y solusi Analitik =');
disp(y2(m));
disp('---------------------------------');
disp('Bola jatuh solusi Numerik =');
disp(yjatuh(j));
disp('Bola jatuh solusi Analitik =');
disp(yjatuh2(n));