%%Read in Datasets - to become a method for trialling different datasets
%Author: Mel Beckerleg
%see also rank_k_problem // rank_k_solve
% input: prob - a problem structure with the following fields
%               - dataset_name 
%                   from {'chembl','jester','RC','MovieLens','Book-Crossings'}
%               -m,n the size of the dataset
%
function [db,prob]=load_dataset(prob)
homepath=prob.homepath;
dataset_name=prob.dataset_name;
m=prob.m;n=prob.n;
seed=prob.random_subset_trial_idx;
%% Load matrix
%format: output a M by N sparse {+1,-1} matrix with zeros for missing entries. 
if strcmp(dataset_name,'chembl')
    load([homepath 'MP2/Data/read_ch/mats/ch_dense_set'],'supported_ch'); 
    db=supported_ch;
else
    if strcmp(dataset_name,'jester')
        db=readtable([homepath 'Datasets/jester-data-1/jester-data-1.xls']);
        db=table2array(db);
        db=db(2:end,:);
        mat_unknown=db>11;
        db=(db>0)*1-(db<0)*1.;
        db(mat_unknown)=0;
        db=sparse(db);
    else
            if strcmp(dataset_name,'RC')
                db=readtable([homepath 'Datasets/RCdata/rating_final.csv']);
                [vals,idx1,idx2]=unique(db(:,1));num=size(vals);num=num(1);updated_vals=linspace(1,num,num);
                rows=updated_vals(idx2);
                [vals,idx1,idx2]=unique(db(:,2));num=size(vals);num=num(1);updated_vals=linspace(1,num,num);
                cols=updated_vals(idx2);
                vals=table2array(db(:,3))+1;
                vals=2*(vals>2)-1;
                db= sparse(rows,cols,vals,max(rows),max(cols),length(vals));
            else
                if strcmp(dataset_name,'movielens')                    
                    db=table2array(readtable([homepath 'Datasets/ml-100k/udata']));
                    vals1=-(db(:,3)<=3);
                    vals=vals1+(db(:,3)>3);
                    %vals=2*(db(:,3)>3)-1;
                    db=sparse(db(:,1),db(:,2),vals,max(db(:,1)),max(db(:,2)),length(db(:,1)));

                else
                        if strcmp(dataset_name,'BX') 
                            %THIS IS SUUUUPER SPARSE!
                            db=readtable([homepath 'Datasets/BX-CSV_Dump/BX-Book-Ratings.csv']);
                            [vals,idx1,idx2]=unique(db(:,1));num=size(vals);num=num(1);updated_vals=linspace(1,num,num);
                            rows=updated_vals(idx2);
                            [vals,idx1,idx2]=unique(db(:,2));num=size(vals);num=num(1);updated_vals=linspace(1,num,num);
                            cols=updated_vals(idx2);
                            vals=str2double(table2array(db(:,3)));
                            vals=2*(vals>5)-1;
                            db= sparse(rows,cols,vals,max(rows),max(cols),length(vals));
                            db=db(find((sum(db')>5)),:);db=db(:,find((sum(db)>5)));db=db(find((sum(db')>5)),:);
                            size(db)
                        else
                            if strcmp(dataset_name,'reddit')
                                
                            rows=csvread('/home/user/Documents/Mel/Ethera/Datasets/Reddit/reddit_more_than_10_rows.csv');
                            cols=csvread('/home/user/Documents/Mel/Ethera/Datasets/Reddit/reddit_more_than_10_cols.csv');
                            vals=csvread('/home/user/Documents/Mel/Ethera/Datasets/Reddit/reddit_more_than_10_vals.csv');
                            db=sparse(rows+1,cols+1,vals);

                            else
                                if strcmp(dataset_name,'netflix')
                                   db=csvread('/home/user/Documents/Mel/Ethera/Datasets/Netflix/ml_csv.csv',1,1);
                                   db=sparse(db);
                                   db(find(db==1))=-1;db(find(db==2))=-1;db(find(db==3))=0;
                                   db(db>3)=1;
                                else
                                disp('Unknown dataset')
                            
                                end
                            end
                        end
                        
                end
            end
    end
        
end

%% select a smaller subset (& drop zero rows/columns)
[M,N]=size(db);
m_val=min(M,m);n_val=min(N,n);
rng(seed);
row_idx=randi(M,[m_val,1]);
col_idx=randi(N,[n_val,1]);
db=db(row_idx,:);db=db(:,col_idx);

% Remove zero rows
db( all(~db,2), : ) = [];
% Remove zero columns
db( :, all(~db,1) ) = [];
% Remove zero rows
db( all(~db,2), : ) = [];

[prob.m, prob.n]=size(db);   
end


    
   

