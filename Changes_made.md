Log book in code changes and modifications
=============================
---

This notebook is a history of the changes made to the code to keep track of these.

## Compilation issue on AMD systems Mi300 and Lumi.
### 06 Apr 25 : 12:20

The compilation issue is corrected with the following modification:

```
[Linux][10:59:03] dc-bonn2@a89612a82902:~/SwanSea/SourceCodes/Grid-DWF-Telos/Grid (develop)() =>$ git diff  Grid/threads/Accelerator.h
diff --git a/Grid/threads/Accelerator.h b/Grid/threads/Accelerator.h
index b2a40e7b..2c390d66 100644
--- a/Grid/threads/Accelerator.h
+++ b/Grid/threads/Accelerator.h
@@ -518,7 +518,8 @@ inline void acceleratorMemSet(void *base,int value,size_t bytes) { auto discard=

inline void acceleratorCopyDeviceToDeviceAsynch(const void *from,void *to,size_t bytes) // Asynch
{
-  auto discard=hipMemcpyDtoDAsync(to,from,bytes, copyStream);
+       void* from_c = const_cast<void*>(from);
+  auto discard=hipMemcpyDtoDAsync(to,from_c,bytes, copyStream);
   }
   inline void acceleratorCopyToDeviceAsync(const void *from, void *to, size_t bytes, hipStream_t stream = copyStream) {
   auto r = hipMemcpyAsync(to,from,bytes, hipMemcpyHostToDevice, stream);
[Linux][10:59:16] dc-bonn2@a89612a82902:~/SwanSea/SourceCodes/Grid-DWF-Telos/Grid (develop)() =>$
```


## Compilation issue on AMD systems Mi300 and Lumi.
### 06 Apr 25 : 12:20

update of the hipMallocHost 


```
inline void *acceleratorAllocHost(size_t bytes)
{
  void *ptr=NULL;
  auto err = hipMallocHost((void **)&ptr,bytes);
  //auto err = hipHostMalloc((void **)&ptr,bytes);

  if( err != hipSuccess ) {
    ptr = (void *) NULL;
    fprintf(stderr," hipMallocManaged failed for %ld %s \n",bytes,hipGetErrorString(err)); fflush(stderr);
  }
  return ptr;
};
```
