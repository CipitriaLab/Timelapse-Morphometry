function digit = get_digit_3(number)
% By default this is 4 digits!
if number == 0
    digit = strcat('00',num2str(number));
elseif number < 10
    digit = strcat('00',num2str(number));
elseif number <100
    digit = strcat('0',num2str(number));
elseif number < 999
        digit = num2str(number);
else
    error('The number is too high for 3 digits');
%     digit = num2str(number);
end
end