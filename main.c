#include "return_codes.h"

#include <malloc.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>

#define MIN(i, j) (((i) < (j)) ? (i) : (j))
#define ABS(i) (((i) >= 0) ? (i) : (-i))

void get_message(int code)
{
	switch (code)
	{
	case 0:
		fprintf(stderr, "The operation completed successfully");
		break;
	case 1:
		fprintf(stderr, "File can't be opened");
		break;
	case 2:
		fprintf(stderr, "Not enough memory, memory allocation failed");
		break;
	case 3:
		fprintf(stderr, "The data is invalid");
		break;
	case 4:
		fprintf(stderr, "The parameter or number of parameters (argv) is incorrect");
		break;
	case 5:
		fprintf(stderr, "Unsupported functionality");
		break;
	default:
		fprintf(stderr, "Unknown error");
	}
}
int write_matrix_to_file(double *eigenvalues, double *complex_eigenvalues, char *File, int n)
{
	FILE *outputFile = fopen(File, "w");
	if (outputFile == NULL)
		return ERROR_CANNOT_OPEN_FILE;
	for (int i = 0; i < n; ++i)
	{
		if (fprintf(outputFile, "%g", eigenvalues[i]) < 0)
		{
			fclose(outputFile);
			return ERROR_UNKNOWN;
		}
		if (complex_eigenvalues[i] == 0)
		{
			if (fprintf(outputFile, "\n") < 0)
			{
				fclose(outputFile);
				return ERROR_UNKNOWN;
			}
		}
		else if (complex_eigenvalues[i] > 0)
		{
			if (fprintf(outputFile, " +%gi\n", complex_eigenvalues[i]) < 0)
			{
				fclose(outputFile);
				return ERROR_UNKNOWN;
			}
		}
		else
		{
			if (fprintf(outputFile, " %gi\n", complex_eigenvalues[i]) < 0)
			{
				fclose(outputFile);
				return ERROR_UNKNOWN;
			}
		}
	}
	fclose(outputFile);
	return SUCCESS;
}

int read_matrix_from_file(char *File, int *n, double **point_on_matrix)
{
	FILE *inputFile = fopen(File, "r");
	if (inputFile == NULL)
	{
		return ERROR_CANNOT_OPEN_FILE;
	}
	if (fscanf(inputFile, "%d", n) != 1)
	{
		fclose(inputFile);
		return ERROR_DATA_INVALID;
	}
	if (*n < 1)
	{
		fclose(inputFile);
		return ERROR_DATA_INVALID;
	}
	double *matrix = malloc(sizeof(double) * *n * *n);
	*point_on_matrix = matrix;
	if (matrix == NULL)
	{
		free(matrix);
		fclose(inputFile);
		return ERROR_OUT_OF_MEMORY;
	}
	for (int i = 0; i < *n; ++i)
	{
		for (int j = 0; j < *n; ++j)
		{
			if (fscanf(inputFile, "%lf", &matrix[i * *n + j]) != 1)
			{
				free(matrix);
				fclose(inputFile);
				return ERROR_DATA_INVALID;
			}
		}
	}
	fclose(inputFile);
	return 0;
}
double norma_Vector(const double *vector, int len)
{
	double sum = 0;
	for (int i = 0; i < len; ++i)
	{
		sum += vector[i] * vector[i];
	}
	return sqrt(sum) != 0 ? sqrt(sum) : 1;
}
double norma_Vector_in_Matrix(const double *matrix, int y, int x, int end, int n)
{
	double sum = 0;
	for (int i = y; i <= end; ++i)
	{
		sum += matrix[i * n + x] * matrix[i * n + x];
	}
	return sqrt(sum) != 0 ? sqrt(sum) : 1;
}
double norma_matrix(const double *matrix, int n)
{
	double sum = 0;
	for (int i = 0; i < n; ++i)
	{
		for (int j = 0; j < n; ++j)
		{
			sum += matrix[i * n + j] * matrix[i * n + j];
		}
	}
	return sqrt(sum);
}
void multiply_V_N(double *vector, double k, int start, int end)
{
	for (int i = start; i <= end; ++i)
	{
		vector[i - start] *= k;
	}
}

void multiply_V_M(const double *uVector1, const double *HessMatrix, int start1, int end1, int start2, int end2, int n, double *uVector2, bool V_first)
{
	for (int i = start2; i <= end2; ++i)
	{
		double sum = 0;
		for (int j = start1; j <= end1; ++j)
		{
			sum += uVector1[j - start1] * HessMatrix[V_first ? (j * n + i) : (i * n + j)];
		}
		uVector2[i - start2] = sum;
	}
}

void subtract_M_M(double *HessMatrix, const double *Householder, int start1, int end1, int start2, int end2, int n)
{
	for (int i = start1; i <= end1; ++i)
	{
		for (int j = start2; j <= end2; ++j)
		{
			HessMatrix[i * n + j] -= Householder[i * n + j];
		}
	}
}
void get_uVector_from_vector(double *uVector, const double *Vector, int len)
{
	double p = Vector[0] >= 0 ? 1 : -1;
	double norma_X = norma_Vector(Vector, len);
	uVector[0] = Vector[0] + p * norma_X;
	for (int i = 1; i < len; ++i)
	{
		uVector[i] = Vector[i];
	}
	double norma_U = norma_Vector(uVector, len);
	for (int i = 0; i < len; ++i)
	{
		uVector[i] /= norma_U;
	}
}

void get_uVector_from_matrix(double *uVector, const double *HessMatrix, int start, int end, int n)
{
	double p = HessMatrix[start * n + start - 1] >= 0 ? 1 : -1;
	double norma_X = norma_Vector_in_Matrix(HessMatrix, start, start - 1, end, n);
	uVector[0] = HessMatrix[start * n + start - 1] + p * norma_X;
	for (int i = start + 1; i <= end; ++i)
	{
		uVector[i - start] = HessMatrix[i * n + start - 1];
	}
	double norma_U = norma_Vector(uVector, end - start + 1);
	for (int i = start; i <= end; ++i)
	{
		uVector[i - start] /= norma_U;
	}
}

void VxV(const double *uVector1, const double *uVector2, double *Householder, int start1, int end1, int start2, int end2, int n)
{
	for (int i = start1; i <= end1; ++i)
	{
		for (int j = start2; j <= end2; ++j)
		{
			Householder[i * n + j] = uVector1[i - start1] * uVector2[j - start2];
		}
	}
}

void Householder_reflector(
	double *uVector1,
	double *uVector2,
	double *Householder,
	double *HessMatrix,
	int fi_start1,
	int fi_end1,
	int fi_start2,
	int fi_end2,
	int se_start1,
	int se_end1,
	int se_start2,
	int se_end2,
	int n)
{
	multiply_V_M(uVector1, HessMatrix, fi_start1, fi_end1, fi_start2, fi_end2, n, uVector2, true);
	multiply_V_N(uVector1, 2, 0, fi_end1 - fi_start1);
	VxV(uVector1, uVector2, Householder, fi_start1, fi_end1, fi_start2, fi_end2, n);
	subtract_M_M(HessMatrix, Householder, fi_start1, fi_end1, fi_start2, fi_end2, n);

	for (int i = 0; i <= fi_end1 - fi_start1; ++i)
		uVector1[i] /= 2;
	multiply_V_M(uVector1, HessMatrix, se_start2, se_end2, se_start1, se_end1, n, uVector2, false);
	multiply_V_N(uVector1, 2, 0, fi_end1 - fi_start1);
	VxV(uVector2, uVector1, Householder, se_start1, se_end1, se_start2, se_end2, n);
	subtract_M_M(HessMatrix, Householder, se_start1, se_end1, se_start2, se_end2, n);
}

int get_Hessenberg_form(double *HessMatrix, int n)
{
	double *uVector1 = malloc(sizeof(double) * n);
	double *uVector2 = malloc(sizeof(double) * n);
	double *Householder = malloc(sizeof(double) * n * n);
	if (uVector1 == NULL || uVector2 == NULL || Householder == NULL)
	{
		free(uVector1);
		free(uVector2);
		free(Householder);
		return ERROR_OUT_OF_MEMORY;
	}
	for (int i = 0; i < n - 2; ++i)
	{
		get_uVector_from_matrix(uVector1, HessMatrix, i + 1, n - 1, n);
		Householder_reflector(uVector1, uVector2, Householder, HessMatrix, i + 1, n - 1, i, n - 1, 0, n - 1, i + 1, n - 1, n);
	}
	free(uVector1);
	free(uVector2);
	free(Householder);

	return SUCCESS;
}

void get_b_c(const double *HessMatrix, int size, int n, double *b, double *c)
{
	double trail = HessMatrix[(size - 2) * n + size - 2] + HessMatrix[(size - 1) * n + size - 1];
	double determinant =
		HessMatrix[(size - 2) * n + size - 2] * HessMatrix[(size - 1) * n + size - 1] -
		HessMatrix[(size - 2) * n + size - 1] * HessMatrix[(size - 1) * n + size - 2];
	if (trail * trail > 4 * determinant)
	{
		double first_x = (trail + sqrt(trail * trail - 4 * determinant)) / 2;
		double second_x = (trail - sqrt(trail * trail - 4 * determinant)) / 2;
		if (ABS(first_x - HessMatrix[(size - 1) * n + size - 1]) < ABS(second_x - HessMatrix[(size - 1) * n + size - 1]))
			second_x = first_x;
		else
			first_x = second_x;
		*b = -first_x - second_x;
		*c = first_x * second_x;
	}
	else
	{
		*b = -trail;
		*c = determinant;
	}
}

int QR_Francis_double_step(double *HessMatrix, int n)
{
	double *uVector1 = malloc(sizeof(double) * n);
	double *uVector2 = malloc(sizeof(double) * n);
	double *Householder = malloc(sizeof(double) * n * n);
	if (uVector1 == NULL || uVector2 == NULL || Householder == NULL)
	{
		free(uVector1);
		free(uVector2);
		free(Householder);
		return ERROR_OUT_OF_MEMORY;
	}

	int cur_size = n;
	double min_val = norma_matrix(HessMatrix, n) * 1e-8;
	while (cur_size > 2)
	{
		if (ABS(HessMatrix[n * (cur_size - 1) + cur_size - 2]) < min_val)
		{
			HessMatrix[n * (cur_size - 1) + cur_size - 2] = 0;
			cur_size--;
		}
		else if (ABS(HessMatrix[n * (cur_size - 2) + cur_size - 3]) < min_val)
		{
			HessMatrix[n * (cur_size - 2) + cur_size - 3] = 0;
			cur_size--;
		}
		else
		{
			double b;
			double c;
			get_b_c(HessMatrix, cur_size, n, &b, &c);
			uVector2[0] = HessMatrix[0] * HessMatrix[0] + HessMatrix[1] * HessMatrix[n] + b * HessMatrix[0] + c;
			uVector2[1] = HessMatrix[n] * HessMatrix[0] + HessMatrix[n + 1] * HessMatrix[n] + b * HessMatrix[n];
			uVector2[2] = HessMatrix[2 * n] * HessMatrix[0] + HessMatrix[2 * n + 1] * HessMatrix[n];
			get_uVector_from_vector(uVector1, uVector2, 3);
			Householder_reflector(uVector1, uVector2, Householder, HessMatrix, 0, 2, 0, cur_size - 1, 0, cur_size - 1, 0, 2, n);
			for (int i = 0; i <= cur_size - 3; ++i)
			{
				int r = MIN(i + 3, cur_size - 1);
				get_uVector_from_matrix(uVector1, HessMatrix, i + 1, r, n);
				Householder_reflector(uVector1, uVector2, Householder, HessMatrix, i + 1, r, 0, cur_size - 1, 0, cur_size - 1, i + 1, r, n);
				HessMatrix[r * n + i] = 0;
			}
		}
	}
	free(uVector1);
	free(uVector2);
	free(Householder);
	return SUCCESS;
}

void get_eigenvalues(double *HessMatrix, double *eigenvalues, double *complex_eigenvalues, int n)
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if (ABS(HessMatrix[i * n + j]) < 1e-10)
			{
				HessMatrix[i * n + j] = 0;
			}
		}
	}
	int count = 0;
	while (count < n)
	{
		if (count == n - 1 || HessMatrix[n * (count + 1) + count] == 0)
		{
			eigenvalues[count] = HessMatrix[count * n + count];
			complex_eigenvalues[count] = 0;
			count++;
		}
		else
		{
			double b = -(HessMatrix[count * n + count] + HessMatrix[(count + 1) * n + count + 1]);
			double c = HessMatrix[count * n + count] * HessMatrix[(count + 1) * n + count + 1] -
					   HessMatrix[count * n + count + 1] * HessMatrix[(count + 1) * n + count];
			double discriminant = b * b - 4 * c;
			if (discriminant >= 0)
			{
				if (HessMatrix[count * n + count] >= 0)
				{
					eigenvalues[count] = (-b + sqrt(discriminant)) / 2;
					eigenvalues[count + 1] = (-b - sqrt(discriminant)) / 2;
				}
				else
				{
					eigenvalues[count] = (-b - sqrt(discriminant)) / 2;
					eigenvalues[count + 1] = (-b + sqrt(discriminant)) / 2;
				}
				complex_eigenvalues[count] = 0;
				complex_eigenvalues[count + 1] = 0;
			}
			else
			{
				eigenvalues[count] = -b / 2;
				eigenvalues[count + 1] = -b / 2;
				complex_eigenvalues[count] = sqrt(-discriminant) / 2;
				complex_eigenvalues[count + 1] = -sqrt(-discriminant) / 2;
			}
			count += 2;
		}
	}
}

int main(int argc, char *argv[])
{
	if (argc != 3)
	{
		fprintf(stderr, "The parameter or number of parameters (argv) is incorrect\nfirst argument: input file\nsecond "
			   "argument: output file\n");
		return ERROR_PARAMETER_INVALID;
	}
	int n;
	double **point_on_matrix[1];
	int return_code = read_matrix_from_file(argv[1], &n, (double **)point_on_matrix);
	double *matrix = (double *)*point_on_matrix;
	if (return_code != SUCCESS)
	{
		get_message(return_code);
		return return_code;
	}

	double *HessMatrix = matrix;
	return_code = get_Hessenberg_form(HessMatrix, n);
	if (return_code != SUCCESS)
	{
		free(HessMatrix);
		get_message(return_code);
		return return_code;
	}

	return_code = QR_Francis_double_step(HessMatrix, n);
	if (return_code != SUCCESS)
	{
		free(HessMatrix);
		get_message(return_code);
		return return_code;
	}

	double *eigenvalues = malloc(sizeof(double) * n);
	double *complex_eigenvalues = malloc(sizeof(double) * n);
	if (eigenvalues == NULL || complex_eigenvalues == NULL)
	{
		free(HessMatrix);
		free(eigenvalues);
		free(complex_eigenvalues);
		get_message(ERROR_OUT_OF_MEMORY);
		return ERROR_OUT_OF_MEMORY;
	}
	get_eigenvalues(HessMatrix, eigenvalues, complex_eigenvalues, n);

	return_code = write_matrix_to_file(eigenvalues, complex_eigenvalues, argv[2], n);
	if (return_code != SUCCESS)
	{
		free(HessMatrix);
		free(eigenvalues);
		free(complex_eigenvalues);
		get_message(return_code);
		return return_code;
	}
	return SUCCESS;
}