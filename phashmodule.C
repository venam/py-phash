/*
 * Read: http://www.python.org/dev/peps/pep-0007/
 */

///--- Library needed ---///
#include <Python.h>
#include <structmember.h>
#include <pHash.h>
#include <audiophash.h>

///--- Documentation ---///
/* The module doc string */
static char module_docstring[] =
"Python bindings for libpHash (http://phash.org/)\n\
A perceptual hash is a fingerprint of a multimedia file derived from various \n\
features from its content. Unlike cryptographic hash functions which rely on \n\
the avalanche effect of small changes in input leading to drastic changes in \n\
the output, perceptual hashes are \"close\" to one another if the features are\n\
similar.\n\
";

/* The functions doc string */
PyDoc_STRVAR( compare_images__doc__,
"compare_images(file1,file2,pcc=0.0,sigma=3.5,gamma=1.0,N=180,threshold=0.90) -> int\n\n\
Compare 2 images given the file names  \n\
Keyword arguments: \n\
file1   -- char string of first image file \n\
file2   -- char string of second image file \n\
pcc   -- (out) double value for peak of cross correlation (default 0.0)\n\
sigma   -- double value for deviation of gaussian filter (default 3.5)\n\
gamma   -- double value for gamma correction of images (default 1.0)\n\
N       -- int number for number of angles (default 180)\n\
threshold (default 0.90) \n\
return   -- int 0 (false) for different image, 1 (true) for same images, less \n\
than 0 for error \n\
");
PyDoc_STRVAR( imagehash__doc__,
"imagehash(file) -> ulong64 (aka PyLong)\n\n\
Compute dct robust image hash\n\
Keyword arguments: \n\
file     -- string variable for name of file\n\
return   -- hash of type ulong64 (aka PyLong)\n\
");
PyDoc_STRVAR( mh_imagehash__doc__,
"mh_imagehash(filename, alpha=2.0f, lvl=1.0f) -> uint8_t[]\n\n\
create MH image hash for filename image\n\
Keyword arguments: \n\
filename -- string name of image file\n\
alpha   -- int scale factor for marr wavelet (default=2)\n\
lvl   -- int level of scale factor (default = 1)\n\
return   -- uint8_t array\n\
");
PyDoc_STRVAR( image_digest__doc__,
"image_digest(file, sigma=1.0, gamma=1.0, int N=180) -> pHash.Digest\n\n\
Compute the image digest given the file name.(radial hash)\n\
Keyword arguments: \n\
file     -- string value for file name of input image.\n\
sigma   -- double value for the deviation for gaussian filter\n\
gamma   -- double value for gamma correction on the input image.\n\
N       -- int value for number of angles to consider\n\
return   -- a Digest struct\n\
");
PyDoc_STRVAR( hamming_distance__doc__,
"hamming_distance(hash1,hash2) -> int\n\n\
Compute the hamming distance between two hash (dct)\n\
Keyword arguments: \n\
hash1   -- The first hash (ulong64)\n\
hash2   -- The second hash (ulong64)\n\
return   -- The distance (over 100)\n\
");
PyDoc_STRVAR( hamming_distance2__doc__,
"hamming_distance2(hashA[], hashB[]) -> double\n\n\
Compute hamming distance between two byte arrays (Mexican Hat)\n\
Keyword arguments: \n\
hashA   -- byte array for first hash\n\
hashB   -- byte array for second hash\n\
return   -- double value for normalized hamming distance\n\
");
PyDoc_STRVAR( crosscorr__doc__,
"crosscorr(x, y, threshold=0.90)\n\n\
Compute the cross correlation of two series vectors (image_digest)\n\
Keyword arguments: \n\
x       -- Digest struct\n\
y       -- Digest struct\n\
threshold-- double value for the threshold value for which 2 images are\n\
considered the same or different.\n\
return   -- (ret,pcc)\n\
ret   -- int value - 1 (true) for same, 0 (false) for different, < 0 for error\n\
pcc   -- double value the peak of cross correlation\n\
");
PyDoc_STRVAR( compare_audio__doc__,
"compare_audio(file1, file2, sr=8000, channels=1, threshold=0.30, block_size=256) -> double\n\n\
Compare 2 audio files audio hash distance\n\
Keyword arguments: \n\
sr       -- sample rate on which to base the audiohash\n\
channels -- nb channels to convert to (always 1) unused\n\
threshold-- threshold value to compare successive blocks, 0.25, 0.30, 0.35\n\
block_size--length of block_size, 256\n\
return   -- confidence score\n\
");
PyDoc_STRVAR( compare_texthash__doc__,
"compare_texthash(file1, file2) -> pHash.TxtCompare \n\n\
textual hash for file\n\
Keyword arguments: \n\
file1    -- name of the first file\n\
file2    -- name of the second file\n\
return   -- A struct containing the hash length of the first file,  the hash \n\
length of the second file, the number of matches, and the maximum consecutive \n\
similar hashes.\n\
");
///--- Globals ---///
static PyObject *pHashError;

typedef struct pHashDigest {
    PyObject_HEAD
    PyObject *id;              /* hash id */
    PyObject *coeffs;          /* the head of the digest integer coefficient array */
    int size;
} pHashDigest;

typedef struct pHashTextCompare {
    PyObject_HEAD
    int length1;
    int length2;
    int matches;
    int maxlength;
} pHashTextCompare;

static PyTypeObject pHashDigestType = {
    PyObject_HEAD_INIT(NULL)
};

static PyTypeObject pHashTextCompareType = {
    PyObject_HEAD_INIT(NULL)
};

static PyMemberDef pHashDigest_members[] = {
    {(char*)"id", T_STRING, offsetof(pHashDigest, id), 0, (char*)"id"},
    {(char*)"coeffs", T_OBJECT, offsetof(pHashDigest, coeffs), 0, (char*)"coeffs"},
    {(char*)"size", T_INT, offsetof(pHashDigest, size), 0, (char*)"size"},
    {NULL}
};
static PyMemberDef pHashTextCompare_members[] = {
    {(char*)"length1", T_INT, offsetof(pHashTextCompare, length1), 0, (char*)"length1"},
    {(char*)"length2", T_INT, offsetof(pHashTextCompare, length2), 0, (char*)"length2"},
    {(char*)"matches", T_INT, offsetof(pHashTextCompare, matches), 0, (char*)"matches"},
    {(char*)"maxlength", T_INT, offsetof(pHashTextCompare, maxlength), 0, (char*)"maxlength"},
    {NULL}
};

///--- Foo Prototypes ---///
static PyObject * PyList_FromUint8Array(uint8_t *array, int len);
static uint8_t  * arrayFromPyList(PyObject* pyList);
static bool file_ready_for_reading (const char *filename);

static PyObject * phash_compare_images(PyObject *self, PyObject *args, PyObject *keywds);
static PyObject * phash_imagehash(PyObject *self, PyObject *args);
static PyObject * phash_mh_imagehash(PyObject *self, PyObject *args, PyObject *keywds);
static PyObject * phash_image_digest(PyObject *self, PyObject *args, PyObject *keywds);
static PyObject * phash_hamming_distance(PyObject *self, PyObject *args);
static PyObject * phash_hamming_distance2(PyObject *self, PyObject *args);
static PyObject * phash_crosscorr(PyObject *self, PyObject *args, PyObject *keywds);
static PyObject * phash_compare_texthash(PyObject *self, PyObject *args);
static PyObject * phash_compare_audio(PyObject *self, PyObject *args, PyObject *keywds);


///--- Definition of the pHash pythons library methods ---///
static PyMethodDef pHash_methods[] = {
    { "compare_images", (PyCFunction)phash_compare_images, METH_VARARGS|METH_KEYWORDS
        , compare_images__doc__ },
    { "imagehash", phash_imagehash, METH_VARARGS,
        imagehash__doc__ },
    { "mh_imagehash", (PyCFunction)phash_mh_imagehash, METH_VARARGS|METH_KEYWORDS,
        mh_imagehash__doc__},
    { "image_digest", (PyCFunction)phash_image_digest, METH_VARARGS|METH_KEYWORDS,
        image_digest__doc__ },
    { "hamming_distance", phash_hamming_distance, METH_VARARGS,
        hamming_distance__doc__ },
    { "hamming_distance2", phash_hamming_distance2, METH_VARARGS,
        hamming_distance2__doc__ },
    { "crosscorr", (PyCFunction)phash_crosscorr, METH_VARARGS|METH_KEYWORDS,
        crosscorr__doc__ },
    { "compare_texthash", phash_compare_texthash, METH_VARARGS,
        compare_texthash__doc__ },
    { "compare_audio", (PyCFunction)phash_compare_audio, METH_VARARGS|METH_KEYWORDS,
        compare_audio__doc__ },
    { NULL, NULL, 0, NULL}
};

///--- Useful functions for the module (not in it) --///

/* Return a new PyList from a byte array (uint8_t / unsigned char array) */
static PyObject *
PyList_FromUint8Array(uint8_t *array, int len)
{
    PyObject *pyList;
    int i;
    if (len < 0) return NULL;

    pyList = PyList_New(len);
    for (i = 0; i < len; i++) {
        PyObject *item = PyInt_FromLong((long) array[i]);
        PyList_SetItem(pyList, i, item);
    }
    return pyList;
}

/* Returns a byte (uint8_t) array from a PyList */
static uint8_t*
arrayFromPyList(PyObject* pyList)
{
    uint8_t *newarr = NULL;
    Py_ssize_t i;
    Py_ssize_t len;

    if (!pyList) return NULL;

    len = PyList_Size(pyList);
    if (len < 0) return NULL;

    newarr = (uint8_t*) malloc(len * sizeof(uint8_t));
    if (!newarr) return NULL;

    for (i = 0; i < len; i++)
        newarr[i] = (uint8_t) PyInt_AsLong(PyList_GetItem(pyList, i));

    return newarr;
}

/* return false if the file doesn't exist or does not have permission to read */
static bool
file_ready_for_reading (const char *filename)
{
    if( access( filename, R_OK ) != 0 ) {
        PyErr_SetString(pHashError,
            "The File you specified does not exist or cannot be read.");
        return false;
    }
    return true;
}

///--- Module functions ---///
/* initialization */
PyMODINIT_FUNC
initpHash(void)
{
    PyObject *m = Py_InitModule3("pHash", pHash_methods, module_docstring);
    if (m ==NULL) return;
    /* Error handler */
    pHashError = PyErr_NewException((char*)"pHash.error", NULL, NULL);
    Py_INCREF(pHashError);
    PyModule_AddObject(m, "error", pHashError);
    /* Digest type */
    pHashDigestType.tp_name           = "pHash.Digest";
    pHashDigestType.tp_basicsize      = sizeof(pHashDigest);
    pHashDigestType.tp_new            = PyType_GenericNew;
    pHashDigestType.tp_methods        = NULL;
    pHashDigestType.tp_members        = pHashDigest_members;
    pHashDigestType.tp_flags          = Py_TPFLAGS_DEFAULT;
    pHashDigestType.tp_doc            = "A pHash radial digest object";
    PyType_Ready(&pHashDigestType);
    Py_INCREF(&pHashDigestType);
    PyModule_AddObject(m, "Digest", (PyObject *)&pHashDigestType);
    /* TextCompare type */
    pHashTextCompareType.tp_name      = "pHash.TextCompare";
    pHashTextCompareType.tp_basicsize = sizeof(pHashTextCompare);
    pHashTextCompareType.tp_new       = PyType_GenericNew;
    pHashTextCompareType.tp_methods   = NULL;
    pHashTextCompareType.tp_members   = pHashTextCompare_members;
    pHashTextCompareType.tp_flags     = Py_TPFLAGS_DEFAULT;
    pHashTextCompareType.tp_doc       = "A structure containing the result of textual hash";
    PyType_Ready(&pHashTextCompareType);
    Py_INCREF(&pHashTextCompareType);
    PyModule_AddObject(m, "TextCompare", (PyObject *)&pHashTextCompareType);
}

static PyObject *
phash_compare_images(PyObject *self, PyObject *args, PyObject *keywds)
{
    /* set keywords and default args */
    static char *kwlist[] = {
        (char*)"file1", (char*)"file2", (char*)"pcc", (char*)"sigma",
        (char*)"gamma", (char*)"N", (char*)"threshold", (char*)NULL};
    const char *file1;
    const char *file2;
    double pcc     = 0.0;
    double sigma     = 3.5;
    double gamma     = 1.0;
    int N           = 180;
    double threshold = 0.90;
    /* take the args from the user */
    if(!PyArg_ParseTupleAndKeywords(args, keywds , "ss|dddid", kwlist,
        &file1, &file2, &pcc, &sigma, &gamma, &N, &threshold))
        return NULL;
    /* test if the files exists and ready for reading */
    if (!file_ready_for_reading(file1)) return NULL;
    if (!file_ready_for_reading(file2)) return NULL;
    /* do the comparison and return the value to the user */
    ph_compare_images(file1, file2, pcc, sigma, gamma, N, threshold);
    return PyFloat_FromDouble(pcc);
}

static PyObject *
phash_imagehash(PyObject *self, PyObject *args)
{
    const char *filename;
    ulong64 hash = 0;
    int test_imagehash;
    if(!PyArg_ParseTuple(args, "s", &filename))
        return NULL;
    /* Check if the file exist and ready for reading */
    if (!file_ready_for_reading(filename)) return NULL;

    test_imagehash = ph_dct_imagehash(filename, hash);
    /* the dct_imagehash has failed */
    if (test_imagehash == -1) {
        PyErr_SetString(pHashError,
            "The DCT robust image hash has failed.");
        return NULL;
    }
    return PyLong_FromUnsignedLongLong(hash);
}

static PyObject *
phash_mh_imagehash(PyObject *self, PyObject *args, PyObject *keywds)
{
    /* set keywords and default args */
    static char *kwlist[] = {
        (char*)"filename", (char*)"alpha", (char*)"lvl", (char*)NULL};
    const char *filename;
    int N        = 0;
    float alpha   = 2.0f;
    float lvl    = 1.0f;
    uint8_t* hash = 0;
    /* take the args from the user */
    if(!PyArg_ParseTupleAndKeywords(args, keywds , "s|ff", kwlist,
        &filename, &alpha, &lvl))
        return NULL;
    /* Check if the file exist and ready for reading */
    if (!file_ready_for_reading(filename)) return NULL;

    hash = ph_mh_imagehash(filename, N, alpha, lvl);
    return PyList_FromUint8Array(hash, N);
}

static PyObject *
phash_image_digest(PyObject *self, PyObject *args, PyObject *keywds)
{
    /* set keywords and default args */
    static char *kwlist[] = {
        (char*)"file", (char*)"sigma", (char*)"gamma", (char*)"N", (char*)NULL};
    const char *filename;
    double sigma=1.0, gamma=1.0;
    Digest dig;
    int N = 180, i;
    PyObject *coeffs, *coeff;
    pHashDigest *phdig;

    if(!PyArg_ParseTupleAndKeywords(args, keywds,"s|ddi:", kwlist,
        &filename, &sigma, &gamma, &N))
        return NULL;
    /* Check if the file exist and ready for reading */
    if (!file_ready_for_reading(filename)) return NULL;

    i = ph_image_digest(filename, sigma, gamma, dig, N);
    if (i<0) {
        PyErr_SetString(pHashError,
            "Computing the image digest of the given file has failed.");
        return NULL;
    }
    phdig = (pHashDigest *)PyObject_New(pHashDigest, &pHashDigestType);

    phdig->id = NULL; // pHash does not make use of this field for radial hashes
    coeffs = PyTuple_New(dig.size);
    for (i=0; i<dig.size; i++) {
        PyTuple_SetItem(coeffs, i, Py_BuildValue("H", dig.coeffs[i]));
    }
    phdig->coeffs = coeffs;
    phdig->size   = dig.size;
    return (PyObject *)phdig;
}

static PyObject *
phash_hamming_distance(PyObject *self, PyObject *args)
{
    ulong64 hash1, hash2;
    int ret;

    if(!PyArg_ParseTuple(args, "KK", &hash1, &hash2))
        return NULL;
    ret = ph_hamming_distance(hash1, hash2);
    return Py_BuildValue("i", ret);
}

static PyObject *
phash_hamming_distance2(PyObject *self, PyObject *args)
{
    PyObject *pyList1, *pyList2;
    uint8_t *uiarr1, *uiarr2;
    int len1, len2;
    double diff = 0.0;

    if (!PyArg_ParseTuple(args, "OO", &pyList1, &pyList2))
        return NULL;
    len1 = PyList_Size(pyList1);
    len2 = PyList_Size(pyList2);

    uiarr1 = arrayFromPyList(pyList1);
    if (!uiarr1) return NULL;
    uiarr2 = arrayFromPyList(pyList2);
    if (!uiarr2) return NULL;

    diff = ph_hammingdistance2(uiarr1, len1, uiarr2, len2);

    free(uiarr1);
    free(uiarr2);
    return PyFloat_FromDouble(diff);
}

static PyObject *
phash_crosscorr(PyObject *self, PyObject *args, PyObject *keywds)
{
    /* set keywords and default args */
    static char *kwlist[] = {
        (char*)"x", (char*)"y", (char*)"threshold", (char*)NULL};
    int ret, size, i;
    double pcc, threshold = 0.90;
    PyObject *py_Digest1, *py_Digest2;
    Digest digest1, digest2;
    uint8_t *coeffs1, *coeffs2;

    if(!PyArg_ParseTupleAndKeywords(args, keywds,"O!O!|d:", kwlist,
        &pHashDigestType, &py_Digest1, &pHashDigestType, &py_Digest2, &threshold))
        return NULL;

    digest1.id = NULL;
    digest2.id = NULL;

    coeffs1 = (uint8_t *)malloc(sizeof(coeffs1) * ((pHashDigest *)py_Digest1)->size);
    coeffs2 = (uint8_t *)malloc(sizeof(coeffs2) * ((pHashDigest *)py_Digest2)->size);
    if(!coeffs1 || !coeffs2){
        PyErr_SetString(pHashError,
            "One of the Digest is malformed.");
        return NULL;
    }

    ret = PyTuple_CheckExact(((pHashDigest *)py_Digest1)->coeffs);
    if(!ret) {
        PyErr_SetString(pHashError,
            "x is not a full Digest.");
        return NULL;
    }
    ret = PyTuple_CheckExact(((pHashDigest *)py_Digest2)->coeffs);
    if(!ret) {
        PyErr_SetString(pHashError,
            "y is not a full Digest.");
        return NULL;
    }

    size = ((pHashDigest *)py_Digest1)->size;
    digest1.size = size;
    digest1.coeffs = coeffs1;
    for (i=0; i<size; i++)
        digest1.coeffs[i] = (uint8_t)PyLong_AsLong(PyTuple_GetItem(((pHashDigest *)py_Digest1)->coeffs, i));

    size = ((pHashDigest *)py_Digest2)->size;
    digest2.size = size;
    digest2.coeffs = coeffs2;
    for (i=0; i<size; i++)
        digest2.coeffs[i] = (uint8_t)PyLong_AsLong(PyTuple_GetItem(((pHashDigest *)py_Digest2)->coeffs, i));

    ret = ph_crosscorr(digest1, digest2, pcc, threshold);

    free(coeffs1);
    free(coeffs2);
    return Py_BuildValue("(i,d)", ret, pcc);
}

static PyObject *
phash_compare_texthash(PyObject *self, PyObject *args)
{

    const char *file1;
    const char *file2;
    TxtHashPoint *hash1, *hash2;
    int nbhashes1 = 0,nbhashes2 = 0, maxlength = 0,count,j=0;
    pHashTextCompare *txtcmp;

    if(!PyArg_ParseTuple(args, "ss", &file1, &file2))
        return NULL;
    /* Check if the file exist and ready for reading */
    if (!file_ready_for_reading(file1) || !file_ready_for_reading(file2))
        return NULL;

    hash1 = ph_texthash(file1,&nbhashes1);
    if (!hash1){
        PyErr_SetString(pHashError,
            "Unable to complete text hash function for file1.");
        return NULL;
    }

    hash2 = ph_texthash(file2,&nbhashes2);
    if (!hash2){
        PyErr_SetString(pHashError,
            "Unable to complete text hash function for file2.");
        return NULL;
    }

    TxtMatch *matches = ph_compare_text_hashes(hash1, nbhashes1, hash2, nbhashes2, &count);
    free(hash1);
    free(hash2);
    if (!matches){
        PyErr_SetString(pHashError,
            "Unable to complete compare function.");
        return NULL;
    }

    /* This is the weirdest thing ever!
     * matches[j].length returns garbage values
     * but the same thing works flawlessly in a separate Cpp file
     * So the texthash works but this value is total garbage.
     */
    for (j=0;j<count;j++){
        if (matches[j].length > maxlength){
            maxlength = matches[j].length;
        }
    }
    free(matches);

    txtcmp = (pHashTextCompare *)PyObject_New(pHashTextCompare, &pHashTextCompareType);
    txtcmp->length1    = nbhashes1;
    txtcmp->length2    = nbhashes2;
    txtcmp->matches    = count;
    txtcmp->maxlength  = maxlength;

    return Py_BuildValue("O", txtcmp);
}

static PyObject *
phash_compare_audio(PyObject *self, PyObject *args, PyObject *keywds)
{
    /* This function crash too much
     * even when running in the example.
     */
    // set keywords and default args
    static char *kwlist[] =
        {(char*)"file1", (char*)"file2", (char*)"sr", (char*)"channels",
        (char*)"threshold", (char*)"block_size", (char*)NULL};

    const char *file1, *file2;
    int sr          = 8000; //sample rate to convert the stream
    int channels    = 1;    //number of channels to convert stream
    float threshold = 0.30; //ber threshold (0.25-0.35)
    int block_size  = 256;  //number of frames to compare at a time
    // take the args from the user
    if(!PyArg_ParseTupleAndKeywords(args, keywds , "ss|iifi", kwlist,
        &file1, &file2, &sr, &channels, &threshold, &block_size))
        return NULL;
    // Check if the file exist and ready for reading
    if (!file_ready_for_reading(file1) || !file_ready_for_reading(file1))
        return NULL;

    float *buf;
    int buflen;
    uint32_t *hash1, *hash2;
    int hashlen1, hashlen2;
    double *cs;
    int Nc;
    int index, i, j;

    buf = ph_readaudio(file1, sr, channels, NULL, buflen);
    if (!buf){
        PyErr_SetString(pHashError,
                "Unable to read audio of file1.");
        return NULL;
    }
    hash1 = ph_audiohash(buf, buflen, sr, hashlen1);
    if (!hash1){
        PyErr_SetString(pHashError,
            "Unable to get hash1.");
        return NULL;
    }
    free(buf);

    buf = ph_readaudio(file2, sr, channels, NULL, buflen);
    if (!buf) {
        PyErr_SetString(pHashError,
            "Unable to read audio of file2.");
        return NULL;
    }
    hash2 = ph_audiohash(buf, buflen, sr, hashlen2);
    if (!hash2) {
        PyErr_SetString(pHashError,
            "Unable to get hash2.");
        return NULL;
    }

    cs = ph_audio_distance_ber(hash1, hashlen1, hash2, hashlen2, threshold, block_size, Nc);
    if (!cs){
        PyErr_SetString(pHashError,
            "Unable to calculate distance.");
        return NULL;
    }

    double max_cs = 0.0;
    for (i=0;i<Nc;i++){
        if (cs[i] > max_cs){
            max_cs = cs[i];
        }
    }

    free(hash2);
    free(hash1);
    free(buf);
    free(cs);
    return Py_BuildValue("d",max_cs);
}
