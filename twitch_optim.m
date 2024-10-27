function twitch_train=twitch_optim(theta,force,Timp,fs)

ST=zeros(size(force));
ST(Timp)=1;

Tlead=0;%round(theta(1));
Tc=theta(1);
Thr=theta(2)*theta(1)/3;
Fmax=theta(3);
Ttw=1000;

[~,F]=RaikovaForceTwitch5p(fs,0,0,Tc,Thr,Ttw,Fmax);

taperWin=tukeywin(Ttw,0.2)';
taperWin(1:Ttw/2)=1;

F=taperWin.*F;

F=[zeros(1,Tlead) F(1:end-Tlead)];

twitch_train=conv(ST,F);
twitch_train((length(force)+1):end)=[];

[b,a] = butter(3,4/(fs/2),'high');
twitch_train=filtfilt(b,a,twitch_train);