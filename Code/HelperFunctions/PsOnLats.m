function L = PsOnLats(ppo,trl,crit)
LVp = ppo(:,2)>0; % positive (p_pre - p_ps)
LVw = ppo(:,3)>crit;
LVpst = [false(trl,1);true(length(ppo)-trl,1)];
L = ppo(LVpst&LVw&LVp,1);