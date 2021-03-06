/*
Fast Artificial Neural Network Library (fann)
Copyright (C) 2003 Steffen Nissen (lukesky@diku.dk)

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include <stdio.h>

#include "fann.h"

int print_callback(unsigned int epochs, float error)
{
	printf("Epochs     %8d. Current MSE-Error: %.10f\n", epochs, error);
	return 0;
}

int main()
{
	const float connection_rate = 1;
	const float learning_rate = (const float)0.7;
	const unsigned int num_layers = 3;
	const unsigned int num_neurons_hidden = 96;
	const float desired_error = (const float)0.001;
	struct fann *ann;
	struct fann_train_data *train_data, *test_data;
	
	unsigned int i = 0;

	printf("Creating network.\n");

	train_data = fann_read_train_from_file("../benchmarks/datasets/robot.train");

	ann = fann_create(connection_rate, learning_rate, num_layers,
		train_data->num_input,
		num_neurons_hidden,
		train_data->num_output);

	printf("Training network.\n");

	fann_set_training_algorithm(ann, FANN_TRAIN_INCREMENTAL);
	
	fann_train_on_data(ann, train_data, 1000, 10, desired_error);
	
	/*fann_train_on_data_callback(ann, data, max_iterations, iterations_between_reports, desired_error, print_callback);*/
	printf("Testing network.\n");

	test_data = fann_read_train_from_file("../benchmarks/datasets/robot.test");

	fann_reset_MSE(ann);
	for(i = 0; i < test_data->num_data; i++){
		fann_test(ann, test_data->input[i], test_data->output[i]);
	}
	printf("MSE error on test data: %f\n", fann_get_MSE(ann));

	printf("Saving network.\n");

	fann_save(ann, "robot_float.net");

	printf("Cleaning up.\n");
	fann_destroy_train(train_data);
	fann_destroy_train(test_data);
	fann_destroy(ann);
	
	return 0;
}
