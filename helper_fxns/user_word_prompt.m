function [ word ] = user_word_prompt( names )
% USER_WORD_PROMPT - is a simple function to get a word from the user
while 1
    word = input('Please enter word: ','s');
    if sum(ismember(names, word)) < 1
        fprintf('OOV\n');
    else
        break;
    end
end

end

