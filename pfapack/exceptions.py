"""PFAPACK: Efficient numerical computation of the Pfaffian."""


class PFAPACKError(Exception):
    """Base exception for PFAPACK errors."""


class InvalidDimensionError(PFAPACKError):
    """Raised when matrix dimensions are invalid."""


class InvalidParameterError(PFAPACKError):
    """Raised when an invalid parameter is provided."""


class ComputationError(PFAPACKError):
    """Raised when the computation fails."""


__all__ = [
    "PFAPACKError",
    "InvalidDimensionError",
    "InvalidParameterError",
    "ComputationError",
]
