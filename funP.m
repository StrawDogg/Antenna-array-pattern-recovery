function res=funP(Phase_Amp,f,dx,dy,HPBW,SLL,Nulls)
    import Config.*
    Phase = Phase_Amp(1,1:Config.N)';
    Amp2 = Phase_Amp(1,Config.N+1:end)';
    [F_dB2] = PlanerArrayPattern(f,Amp2,Phase,dx,dy);
    [HPBW2,SLL2,~,F_dB_Nmlz_Phi0_Ipt] = Get_HPBW_SLL_Nulls(0,F_dB2);
    
    Nulls_cha = Nulls(2,:)-F_dB_Nmlz_Phi0_Ipt(Nulls(1,:));
    nullsCha = sum(Nulls_cha.^2);
    
    % 多SLL指标
    SLL_cha = SLL(2,:) - F_dB_Nmlz_Phi0_Ipt(SLL(1,:));
    SLL_ans = sum(SLL_cha.^2);
    
    % 筛选SLL2,加上有奇效
    shaixuan = find(SLL2 > -20);
    SLL2 = SLL2(shaixuan);
    SLL2 = SLL2 + 20;
    maxSLL2 = max(SLL2);
    sumSLL2 = sum(SLL2) - maxSLL2;
    
    % res = (HPBW-HPBW2)^2+ 0.1 * SLL_ans + 0.5 * nullsCha;
    
    %单独mssa的最佳比例
    res = (HPBW-HPBW2)^2+ 0.3 * SLL_ans + 0.03 * nullsCha + 0.1 * sumSLL2;
    
    % res = (HPBW-HPBW2)^2+ 0.3 * SLL_ans + 0.03 * nullsCha;
    
    % res = (HPBW-HPBW2)^2+(SLL-SLL2)^2 + 0.01 * nullsCha;
    % res = (HPBW-HPBW2)^2+(SLL-SLL2)^2;
%     persistent cnt;
%     if isempty(cnt)
%         cnt=0;
%     else
%         cnt = cnt+1;
%     end
%     fprintf('%d : %f，HPBW:%f,SLL:%f\n',cnt,res,HPBW2,SLL2);
%     disp("!!");
end