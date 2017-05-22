#include "fann.h"

int main(int argc, char* argv[])
{
	struct fann_train_data *train, *test, *all;
	unsigned int i, j;

	if(argc == 2){
		train = fann_read_train_from_file(argv[1]);
		fann_shuffle_train_data(train);
		fann_save_train(train, argv[1]);
		return;
	}
	
	if(argc != 3){
		printf("usage: %s train_file test_file\n", argv[0]);
		printf("will shuffle the data in the two files and save the new data back to them.\n\n");
		printf("usage: %s train_file\n", argv[0]);
		printf("will shuffle the data in the files.\n");
		return -1;
	}

	train = fann_read_train_from_file(argv[1]);
	test = fann_read_train_from_file(argv[2]);
	all = fann_merge_train_data(train, test);

	fann_shuffle_train_data(all);

	for(i = 0; i < train->num_data; i++){
		for(j = 0; j < train->num_input; j++){
			train->input[i][j] = all->input[i][j];
		}
		
		for(j = 0; j < train->num_output; j++){
			train->output[i][j] = all->output[i][j];
		}
	}

	for(i = 0; i < test->num_data; i++){
		for(j = 0; j < test->num_input; j++){
			test->input[i][j] = all->input[i + train->num_data][j];
		}
		
		for(j = 0; j < test->num_output; j++){
			test->output[i][j] = all->output[i + train->num_data][j];
		}
	}

	fann_save_train(train, argv[1]);
	fann_save_train(test, argv[2]);
	
	return 0;
}
