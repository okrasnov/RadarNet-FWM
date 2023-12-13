function visualize_model(Scena)

%========================================================
% v.1.0 - 27.12.2012, OK@TUD
%========================================================



%% =======================================
% Visualization
% =======================================
simb=[{'b*'};{'bo'};{'ro'};{'r*'}]

hf=figure;
maxfig(hf,1);

%subplot(1,2,1)
hold on
for i_radar=1:Scena.N_Radars
    plot(Scena.Radar(i_radar).beam.center_line(:,1),Scena.Radar(i_radar).beam.center_line(:,2),'k--');
    plot([Scena.Radar(i_radar).beam.x_circ,Scena.Radar(i_radar).beam.x_circ(1)],...
        [Scena.Radar(i_radar).beam.y_circ,Scena.Radar(i_radar).beam.y_circ(1)],'k-');
    fill(Scena.Radar(i_radar).beam.x_circ,Scena.Radar(i_radar).beam.y_circ,'r','FaceAlpha',0.1);
end
for ind=1:Scena.N_Highways
    plot(Scena.xlimits,[Scena.Radar(1).position(2)+Scena.Highway(ind).y_Distance,Scena.Radar(1).position(2)+Scena.Highway(ind).y_Distance],'b-');
    if (Scena.Highway(ind).x_min<Scena.Highway(ind).x_max)
        plot([Scena.Highway(ind).x_min,Scena.Highway(ind).x_max],[Scena.Highway(ind).y,Scena.Highway(ind).y],'b-','Linewidth',2);
        for i=1:Scena.Highway(ind).N_Targets
            plot(Scena.Highway(ind).Target(i).position(1),Scena.Highway(ind).Target(i).position(2),char(simb(ind)),'Linewidth',2);
        end
    end
end
for i_radar=1:Scena.N_Radars
    plot(Scena.Radar(i_radar).position(1),Scena.Radar(i_radar).position(2),'r*')
    plot(Scena.ufo.x,Scena.ufo.y,'r.');
    text(Scena.Radar(i_radar).position(1)+10,Scena.Radar(i_radar).position(2)+10,num2str(i_radar),'FontSize',18)
end
hold off
grid on
set(gca,'xlim',Scena.xlimits);
set(gca,'ylim',Scena.ylimits);
axis('square');
set(gca,'fontsize',12);
xlabel('Zonal Distance, relative units','fontsize',14);
ylabel('Meridional Distance, relative units','fontsize',14);

if (Scena.flag_RD_plane==1)
    for i_radar=1:Scena.N_Radars
        i_radar
        for i_Time=1:Scena.N_Times
            i_Time
            hf1=figure;
            %maxfig(hf1,1)
            
            subplot(1,2,1)
            hold on
            for ind_Target=1:Scena.Radar(i_radar).Time(i_Time).N_targets
                plot(Scena.Radar(i_radar).Time(i_Time).Target(ind_Target).Doppler,...
                     Scena.Radar(i_radar).Time(i_Time).Target(ind_Target).Range,...
                     char(simb(1)),'Linewidth',2);
            end
            hold off;
            grid on;
            set(gca,'xlim',[-Scena.Radar(i_radar).Doppler.Umbiguity,Scena.Radar(i_radar).Doppler.Umbiguity]);
            set(gca,'ylim',[0,Scena.Radar(i_radar).Range.Max_Range]);
            axis('square');
            set(gca,'fontsize',12);
            xlabel('Doppler velocity, m/s','fontsize',14);
            ylabel('Range, relative units','fontsize',14);
            title(['Radar #',num2str(i_radar),' Mesurement Time',num2str(Scena.Radar(i_radar).Time(i_Time).Time),' sec'])
            
            subplot(1,2,2)
            imagesc([-Scena.Radar(i_radar).Doppler.Umbiguity:2*Scena.Radar(i_radar).Doppler.Umbiguity/(Scena.Radar(i_radar).Doppler.Bursts):Scena.Radar(i_radar).Doppler.Umbiguity],...
                [1:Scena.Radar(i_radar).Range.Output_Range_Samples].*Scena.Radar(i_radar).Range.Resolution,...
                db(abs(Scena.Radar(i_radar).Time(i_Time).RD_plane)))
            set(gca,'YDir','normal');
            set(gca,'clim',[-20,+30]);
            colorbar
            title(['Radar #',num2str(i_radar),' Mesurement Time',num2str(Scena.Radar(i_radar).Time(i_Time).Time),' sec'])
            disp('Press any key...')
            pause;
            close(hf1)
        end;
    end;
end