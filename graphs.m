clc;  clear all;
f_Lambda=fopen('D:/Praktika/spline/data_Lambda','r');%необходимо добавить адрес файла
m=fgetl(f_Lambda);
i=1;
while ~(feof(f_Lambda))
    m=fgetl(f_Lambda)
    if ~isempty(m)
        m=strrep(m,',','.')
        a(i,:)=str2num(m);
        i=i+1;
    end
end
fclose(f_Lambda);
f_Phi=fopen('D:/Praktika/spline/data_Phi','r');%необходимо добавить адрес файла
m=fgetl(f_Phi);
i=1;
while ~(feof(f_Phi))
    m=fgetl(f_Phi)
    if ~isempty(m)
        m=strrep(m,',','.')
        phi(i,:)=str2num(m);
        i=i+1;
    end
end
fclose(f_Phi);
f_Height=fopen('D:/Praktika/spline/data_Height','r');% необходимо добавить адрес файла
m=fgetl(f_Height);
i=1;
while ~(feof(f_Height))
    m=fgetl(f_Height)
    if ~isempty(m)
        m=strrep(m,',','.')
        h(i,:)=str2num(m);
        i=i+1;
    end
end
fclose(f_Height);
f_Accel=fopen('D:/Praktika/spline/data_Accel','r');% необходимо добавить адрес файла
m=fgetl(f_Accel);
i=1;
while ~(feof(f_Accel))
    m=fgetl(f_Accel)
    if ~isempty(m)
        m=strrep(m,',','.')
        acc(i,:)=str2num(m);
        i=i+1;
    end
end
fclose(f_Accel);
f_Ang=fopen('D:/Praktika/spline/data_Angul','r');% необходимо добавить адрес файла
m=fgetl(f_Ang);
i=1;
while ~(feof(f_Ang))
    m=fgetl(f_Ang)
    if ~isempty(m)
        m=strrep(m,',','.')
        ang(i,:)=str2num(m);
        i=i+1;
    end
end
fclose(f_Ang);

N = i;
i = 1;
for i = 1:N-1
    acc_x = acc(i,1);
    acc_y = acc(i,2);
    acc_z = acc(i,3);
    t = a(i,4);
    %scatter(t, acc_z, 'b');
    %hold on;
end
N = i;
i = 1;
for i = 1:N-1
    ang_x = ang(i,1);
    ang_y = ang(i,2);
    ang_z = ang(i,3);
    t = a(i,4);
    %scatter(t, ang_z, 'b');
    %hold on;
end

N = i;
i = 1;
for i = 1:N
    if(i == 1)
        di = a(i,3)/(6*(a(i,6)-a(i,5)));
        x = meshgrid(a(i,5):0.01:a(i,4));
        y=a(i,1)+a(i,2)*(x-a(i,4)).^1+a(i,3)*(x-a(i,4)).^2/2+di*(x-a(i,4)).^3;
        lambda_itog = y(1,:);
    elseif(i == N)
        di = (a(i-1,3)-a(i-2,3))/(6*(a(i-1,6)-a(i-1,5)));
        x = meshgrid(a(i-1,4):0.01:a(i-1,6));
        y=a(i-1,1)+a(i-1,2)*(x-a(i-1,4)).^1+a(i-1,3)*(x-a(i-1,4)).^2/2+di*(x-a(i-1,4)).^3;
        lambda_itog = [lambda_itog,y(1,:)];
    else
        di = (a(i,3)-a(i-1,3))/(6*(a(i,6)-a(i,5)));
        x = meshgrid(a(i-1,4):0.01:a(i,4));
        y=a(i,1)+a(i,2)*(x-a(i,4)).^1+a(i,3)*(x-a(i,4)).^2/2+di*(x-a(i,4)).^3;
        lambda_itog = [lambda_itog,y(1,:)];
    end
    %plot(x,y*180/3.1415926535,'b');
    %hold on;
end



N = i;
i = 1;
for i = 1:N
    if(i == 1)
        di = phi(i,3)/(6*(phi(i,6)-phi(i,5)));
        x = meshgrid(phi(i,5):0.01:phi(i,4));
        z=phi(i,1)+phi(i,2)*(x-phi(i,4)).^1+phi(i,3)*(x-phi(i,4)).^2/2+di*(x-phi(i,4)).^3;
        phi_itog = z(1,:);
    elseif(i == N)
        di = (phi(i-1,3)-phi(i-2,3))/(6*(phi(i-1,6)-phi(i-1,5)));
        x = meshgrid(phi(i-1,4):0.01:phi(i-1,6));
        z=phi(i-1,1)+phi(i-1,2)*(x-phi(i-1,4)).^1+phi(i-1,3)*(x-phi(i-1,4)).^2/2+di*(x-phi(i-1,4)).^3;
        phi_itog = [phi_itog,z(1,:)];
    else
        di = (phi(i,3)-phi(i-1,3))/(6*(phi(i,6)-phi(i,5)));
        x = meshgrid(phi(i-1,4):0.01:phi(i,4));
        z=phi(i,1)+phi(i,2)*(x-phi(i,4)).^1+phi(i,3)*(x-phi(i,4)).^2/2+di*(x-phi(i,4)).^3;
        phi_itog = [phi_itog,z(1,:)];
    end
    plot(x,z*180/3.1415926535,'b');
    hold on;
end

N = i;
i = 1;
for i = 1:N
    if(i == 1)
        di = h(i,3)/(6*(h(i,6)-h(i,5)));
        x = meshgrid(h(i,5):0.01:h(i,4));
        d=h(i,1)+h(i,2)*(x-h(i,4)).^1+h(i,3)*(x-h(i,4)).^2/2+di*(x-h(i,4)).^3;
        h_itog = d(1,:);
    elseif(i == N)
        di = (h(i-1,3)-h(i-2,3))/(6*(h(i-1,6)-h(i-1,5)));
        x = meshgrid(h(i-1,4):0.01:h(i-1,6));
        d=h(i-1,1)+h(i-1,2)*(x-h(i-1,4)).^1+h(i-1,3)*(x-h(i-1,4)).^2/2+di*(x-h(i-1,4)).^3;
        h_itog = [h_itog,d(1,:)];
    else
        di = (h(i,3)-h(i-1,3))/(6*(h(i,6)-h(i,5)));
        x = meshgrid(h(i-1,4):0.01:h(i,4));
        d=h(i,1)+h(i,2)*(x-h(i,4)).^1+h(i,3)*(x-h(i,4)).^2/2+di*(x-h(i,4)).^3;
        h_itog = [h_itog,d(1,:)];
    end
    %plot(x,d,'b');
    %hold on;
end
%plot(lambda_itog*180/3.1415926535, h_itog);

