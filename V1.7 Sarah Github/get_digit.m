function digit = get_digit(number,varargin)
% digit = get_digit(number,varargin)
% number is the double, varargin  can be:
% varargin = 1 means 1 digit in total
% varargin = 2 means 2 digits in total
% exemple: get_digit(10,3) = '010';
%% This function allows to create a normalized format for naming images.
% there are 4 digits all the times and the first elements is named
if isempty(varargin)
    % then by default it will be 3 digits in total
    digit = get_digit_3(number);
    
else
    nb_digits = varargin{1};
    switch nb_digits
        case 1
            digit = num2str(number);
        case 2
            digit = get_digit_2(number);
        case 3
            digit = get_digit_3(number);
        case 4
            digit = get_digit_4(number);
        otherwise
            error('Only 1 to 4 digits is supported at the moment');
    end
    %     switch nb_digits
    %         case 1
    %             nb_maxi = 9;
    %         case 2
    %             nb_maxi = 99;
    %         case 3
    %             nb_maxi = 999;
    %         case 4
    %             nb_maxi = 9999;
    %         case 5
    %             nb_maxi = 99999;
    %         otherwise
    %             error('Number too high, nto supported yet');
    %     end
    %     switch nb_maxi
    %         case 1
    %             digit = num2str(number);
    %         case 2
    %             digit = strcat('0',num2str(number));
    %         case 3
    %             digit = strcat('00',num2str(number));
    %         case 4
    %             digit = strcat('000',num2str(number));
    %         case 5
    %             digit = strcat('00000',num2str(number));
    %         case 6
    %             digit = strcat('000000',num2str(number));
    %         case 7
    %             digit = strcat('0000000',num2str(number));
    %     end
end



end