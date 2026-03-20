#!/bin/sh

## Generate instantiations of fermion types used in Grid.
## See also README.md

## Implementations to instantiate staggered fermions for
STAG_IMPL_LIST=" \
	   StaggeredImplF \
	   StaggeredImplD  "

## Implementations to instantiate five-dimensional staggered fermions for
STAG5_IMPL_LIST=""

## Implementations to instantiate Wilson fermions for
WILSON_IMPL_LIST=" \
	   WilsonImplF \
	   WilsonImplD \
	   WilsonImplD2 \
	   SpWilsonImplF \
	   SpWilsonImplD \
	   WilsonAdjImplF \
	   WilsonAdjImplD \
	   WilsonTwoIndexSymmetricImplF \
	   WilsonTwoIndexSymmetricImplD \
	   WilsonTwoIndexAntiSymmetricImplF \
	   WilsonTwoIndexAntiSymmetricImplD \
	   SpWilsonTwoIndexAntiSymmetricImplF \
	   SpWilsonTwoIndexAntiSymmetricImplD \
	   SpWilsonTwoIndexSymmetricImplF \
	   SpWilsonTwoIndexSymmetricImplD \
	   GparityWilsonImplF \
	   GparityWilsonImplD "

## Implementations to instantiate compact Wilson fermions for
COMPACT_WILSON_IMPL_LIST=" \
	   WilsonImplF \
	   WilsonImplD "

## Implementations to instantiate domain wall fermions for
DWF_IMPL_LIST=" \
           SpWilsonImplF \
           SpWilsonImplD \
           WilsonImplF \
           WilsonImplD \
           WilsonImplD2 \
           ZWilsonImplF \
           ZWilsonImplD2 "

## Implementations to instantiate G-parity domain wall fermions for
GDWF_IMPL_LIST=" \
	   GparityWilsonImplF \
	   GparityWilsonImplD "

## All implementations to generate
IMPL_LIST="$STAG_IMPL_LIST  $WILSON_IMPL_LIST $DWF_IMPL_LIST $GDWF_IMPL_LIST"

for impl in $IMPL_LIST
do
  echo $impl
  mkdir -p $impl
cat > $impl/impl.h <<EOF
#define IMPLEMENTATION $impl
EOF

done

## Wilson fermion instantiations to generate
CC_LIST="WilsonCloverFermionInstantiation WilsonFermionInstantiation WilsonKernelsInstantiation WilsonTMFermionInstantiation"

for impl in $WILSON_IMPL_LIST
do
for f in $CC_LIST
do
  ln -f -s ../$f.cc.master $impl/$f$impl.cc
done
done

## Compact Wilson clover instantiations to generate
CC_LIST="CompactWilsonCloverFermionInstantiation CompactWilsonCloverFermion5DInstantiation"

for impl in $COMPACT_WILSON_IMPL_LIST
do
for f in $CC_LIST
do
  ln -f -s ../$f.cc.master $impl/$f$impl.cc
done
done

## Domain wall fermion instantiations to generate
CC_LIST=" \
  CayleyFermion5DInstantiation \
  ContinuedFractionFermion5DInstantiation \
  DomainWallEOFAFermionInstantiation  \
  MobiusEOFAFermionInstantiation \
  PartialFractionFermion5DInstantiation \
  WilsonFermion5DInstantiation \
  WilsonKernelsInstantiation "

for impl in $DWF_IMPL_LIST $GDWF_IMPL_LIST
do
for f in $CC_LIST
do
  ln -f -s ../$f.cc.master $impl/$f$impl.cc
done
done

# overwrite the .cc file in Gparity directories
for impl in $GDWF_IMPL_LIST
do
  ln -f -s ../WilsonKernelsInstantiationGparity.cc.master $impl/WilsonKernelsInstantiation$impl.cc
done

## Staggered fermion instantiations to generate
CC_LIST=" \
  ImprovedStaggeredFermion5DInstantiation \
  ImprovedStaggeredFermionInstantiation \
  NaiveStaggeredFermionInstantiation \
  StaggeredKernelsInstantiation "

for impl in $STAG_IMPL_LIST
do
for f in $CC_LIST
do
  ln -f -s ../$f.cc.master $impl/$f$impl.cc
done
done

## Five-dimensional staggered fermion instantiations to **not** generate
CC_LIST=" \
  ImprovedStaggeredFermion5DInstantiation \
  StaggeredKernelsInstantiation "

