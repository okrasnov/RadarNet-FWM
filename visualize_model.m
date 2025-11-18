function visualize_model(Scene)

%========================================================
% v.1.0 - 27.12.2012, OK@TUD
%========================================================



%% =======================================
% Visualization
% =======================================
symb=[{'b*'};{'bo'};{'ro'};{'r*'}]

hf=figure;
maxfig(hf,1);

%subplot(1,2,1)
hold on
for i_radar=1:Scene.N_Radars
    plot(Scene.Radar(i_radar).beam.center_line(:,1),Scene.Radar(i_radar).beam.center_line(:,2),'k--');
    plot([Scene.Radar(i_radar).beam.x_circ,Scene.Radar(i_radar).beam.x_circ(1)],...
        [Scene.Radar(i_radar).beam.y_circ,Scene.Radar(i_radar).beam.y_circ(1)],'k-');
    fill(Scene.Radar(i_radar).beam.x_circ,Scene.Radar(i_radar).beam.y_circ,'r','FaceAlpha',0.1);
end
for ind=1:Scene.N_Highways
    plot(Scene.xlimits,[Scene.Radar(1).position(2)+Scene.Highway(ind).y_Distance,Scene.Radar(1).position(2)+Scene.Highway(ind).y_Distance],'b-');
    if (Scene.Highway(ind).x_min<Scene.Highway(ind).x_max)
        plot([Scene.Highway(ind).x_min,Scene.Highway(ind).x_max],[Scene.Highway(ind).y,Scene.Highway(ind).y],'b-','Linewidth',2);
        for i=1:Scene.Highway(ind).N_Targets
            plot(Scene.Highway(ind).Target(i).position(1),Scene.Highway(ind).Target(i).position(2),char(symb(ind)),'Linewidth',2);
        end
    end
end
for i_radar=1:Scene.N_Radars
    plot(Scene.Radar(i_radar).position(1),Scene.Radar(i_radar).position(2),'r*')
    plot(Scene.ufo.x,Scene.ufo.y,'r.');
    text(Scene.Radar(i_radar).position(1)+10,Scene.Radar(i_radar).position(2)+10,num2str(i_radar),'FontSize',18)
end
hold off
grid on
set(gca,'xlim',Scene.xlimits);
set(gca,'ylim',Scene.ylimits);
axis('square');
set(gca,'fontsize',12);
xlabel('Zonal Distance, relative units','fontsize',14);
ylabel('Meridional Distance, relative units','fontsize',14);

if (Scene.flag_RD_plane==1)
    for i_radar=1:Scene.N_Radars
        i_radar
        for i_Time=1:Scene.N_Times
            i_Time
            hf1=figure;
            %maxfig(hf1,1)
            
            subplot(1,2,1)
            hold on
            for ind_Target=1:Scene.Radar(i_radar).Time(i_Time).N_targets
                plot(Scene.Radar(i_radar).Time(i_Time).Target(ind_Target).Doppler,...
                     Scene.Radar(i_radar).Time(i_Time).Target(ind_Target).Range,...
                     char(symb(1)),'Linewidth',2);
            end
            hold off;
            grid on;
            set(gca,'xlim',[-Scene.Radar(i_radar).Doppler.Umbiguity,Scene.Radar(i_radar).Doppler.Umbiguity]);
            set(gca,'ylim',[0,Scene.Radar(i_radar).Range.Max_Range]);
            axis('square');
            set(gca,'fontsize',12);
            xlabel('Doppler velocity, m/s','fontsize',14);
            ylabel('Range, relative units','fontsize',14);
            title(['Radar #',num2str(i_radar),' Mesurement Time',num2str(Scene.Radar(i_radar).Time(i_Time).Time),' sec'])
            
            subplot(1,2,2)
            imagesc([-Scene.Radar(i_radar).Doppler.Umbiguity:2*Scene.Radar(i_radar).Doppler.Umbiguity/(Scene.Radar(i_radar).Doppler.Bursts):Scene.Radar(i_radar).Doppler.Umbiguity],...
                [1:Scene.Radar(i_radar).Range.Output_Range_Samples].*Scene.Radar(i_radar).Range.Resolution,...
                db(abs(Scene.Radar(i_radar).Time(i_Time).RD_plane)))
            set(gca,'YDir','normal');
            set(gca,'clim',[-20,+30]);
            colorbar
            title(['Radar #',num2str(i_radar),' Mesurement Time',num2str(Scene.Radar(i_radar).Time(i_Time).Time),' sec'])
            disp('Press any key...')
            pause;
            close(hf1)
        end;
    end;
end