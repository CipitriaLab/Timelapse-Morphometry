function digit = get_digit_4(number)
% By default this is 4 digits!
if number == 0
    digit = strcat('000',num2str(number));
elseif number < 10
    digit = strcat('000',num2str(number));
elseif number <100
    digit = strcat('00',num2str(number));
elseif number <1000
    digit = strcat('0',num2str(number));
elseif number < 9999
    digit = num2str(number);
else
    error('The number is too high for 4 digits');
end

end