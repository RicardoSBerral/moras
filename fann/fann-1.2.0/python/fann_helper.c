#include <Python.h>
#include <fann.h>

PyObject *fann_type_to_PyList(fann_type *array,int n)
{
	int i;
	PyObject* res = PyList_New(n);
	for (i = 0; i < n; i++) 
	{
		PyObject *o = PyFloat_FromDouble((double) array[i]);
		PyList_SetItem(res,i,o);
	}

	return res;
}

PyObject *fann_run2(struct fann *ann, fann_type *input)
{
	if (ann == NULL) 
	{
		return NULL;
	}

	return fann_type_to_PyList(fann_run(ann,input),ann->num_output);
}

PyObject *fann_test2(struct fann *ann, fann_type *input, fann_type *desired_output)
{
	if (ann == NULL) 
	{
		return NULL;
	}

	return fann_type_to_PyList(fann_test(ann,input,desired_output),ann->num_output);
}

static PyObject *get_row_from_double_array(fann_type **f, int row, int ncols) 
{
	int i;
	PyObject *result = PyList_New(ncols);

	for (i=0; i<ncols; i++) 
	{
		/* printf("%d:%d %f\n", row, i, f[row][i]); */
		PyList_SetItem(result, i, Py_BuildValue("d", f[row][i]));
	}

	return result;
}

PyObject *get_train_data_input(struct fann_train_data *t, int row) 
{
	if (row < 0 || row >= t->num_data) 
	{
		return Py_BuildValue("");
	}

	return get_row_from_double_array(t->input, row, t->num_input);
}

PyObject *get_train_data_output(struct fann_train_data *t, int row) 
{
	if (row < 0 || row >= t->num_data) 
	{
		return Py_BuildValue("");
	}

	return get_row_from_double_array(t->output, row, t->num_output);
}

