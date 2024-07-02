function s = makeString(v)
% create a comma separated list based on a input vecotr
s = '';
for i = 1:length(v)-1
   s = strcat(s,num2str(v(i)),',');
end
s = strcat(s,num2str(v(end)));
end

