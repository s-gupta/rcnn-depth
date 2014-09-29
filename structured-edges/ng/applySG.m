function szg = applySG(zg, radSG, gtheta)
 	filters = make_filters(radSG, gtheta);
 	for i = 1:length(zg),
 		if(radSG(i) == 0)
 			szg{i} = zg{i};
 		else
 			for o = 1:size(zg{i},3), 
 				szg{i}(:,:,o) = fitparab(abs(zg{i}(:,:,o)),radSG(i),radSG(i)/4,gtheta(o),filters{i,o});
 			end
 		end
 	end
end 
