LZC_table(strcmp(LZC_table.Subject),:).Condition(LZC_table.Subject == 'KET01')




LZC_table.Condition = repmat("Placebo", 762,1);

LZC_table(all(LZC_table.Subject=='KET01',2),:).Condition = deal(temp{:})

LZC_table(all(LZC_table.Subject=='KET01',2) ,:)

ket_lines = boolean(all(LZC_table.Subject=='KET01',2) .* all(LZC_table.Session=='2',2)+...
    all(LZC_table.Subject=='KET02',2) .* all(LZC_table.Session=='1',2)  + ...
    all(LZC_table.Subject=='KET03',2) .* all(LZC_table.Session=='2', 2));


LZC_table(ket_lines,:).Condition =  repmat('Ketamine', 381 ,1);



LZC_completo = [temp_PerConc_merged;LZC_table];

writetable(LZC_table, 'LZC_Baseline_20231018')


