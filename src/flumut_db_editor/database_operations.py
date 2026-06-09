from flumut.flumutdb.models import (
    Annotation,
    Effect,
    Evidence,
    Mapping,
    Marker,
    Mutation,
    Paper,
    Protein,
    Reference,
    Segment,
    Subtype,
)


class ForeignKeyViolationError(Exception):
    """Raised when attempting to delete an item with dependent records."""

    def __init__(self, item_type: str, item_name: str, violations: dict[str, list]):
        self.item_type = item_type
        self.item_name = item_name
        self.violations = violations
        super().__init__(f"Cannot delete {item_type} '{item_name}': {sum(len(v) for v in violations.values())} dependent records")


class DeletionValidator:
    """Validates whether an item can be safely deleted based on FK constraints."""

    @staticmethod
    def can_delete_segment(segment_name: str) -> tuple[bool, dict[str, list]]:
        proteins = list(Protein.select().where(Protein.segment == segment_name))
        references = list(Reference.select().where(Reference.segment == segment_name))
        violations = {}
        if proteins:
            violations['Proteins'] = proteins
        if references:
            violations['References'] = references
        return len(violations) == 0, violations

    @staticmethod
    def can_delete_protein(protein_id: int) -> tuple[bool, dict[str, list]]:
        annotations = list(Annotation.select().where(Annotation.protein_id == protein_id))
        mutations = list(Mutation.select().where(Mutation.protein_id == protein_id))
        violations = {}
        if annotations:
            violations['Annotations'] = annotations
        if mutations:
            violations['Mutations'] = mutations
        return len(violations) == 0, violations

    @staticmethod
    def can_delete_reference(reference_id: int) -> tuple[bool, dict[str, list]]:
        annotations = list(Annotation.select().where(Annotation.reference_id == reference_id))
        mappings = list(Mapping.select().where(Mapping.reference_id == reference_id))
        violations = {}
        if annotations:
            violations['Annotations'] = annotations
        if mappings:
            violations['Mappings'] = mappings
        return len(violations) == 0, violations

    @staticmethod
    def can_delete_mutation(mutation_id: int) -> tuple[bool, dict[str, list]]:
        markers = list(
            Marker.select().join(Marker.mutations.through_model).where(Marker.mutations.through_model.mutation_id == mutation_id)
        )
        mappings = list(Mapping.select().where(Mapping.mutation_id == mutation_id))
        violations = {}
        if markers:
            violations['Markers'] = markers
        if mappings:
            violations['Mappings'] = mappings
        return len(violations) == 0, violations

    @staticmethod
    def can_delete_marker(marker_id: int) -> tuple[bool, dict[str, list]]:
        evidences = list(Evidence.select().where(Evidence.marker_id == marker_id))
        violations = {}
        if evidences:
            violations['Evidences'] = evidences
        return len(violations) == 0, violations

    @staticmethod
    def can_delete_paper(paper_id: int) -> tuple[bool, dict[str, list]]:
        evidences = list(Evidence.select().where(Evidence.paper_id == paper_id))
        violations = {}
        if evidences:
            violations['Evidences'] = evidences
        return len(violations) == 0, violations

    @staticmethod
    def can_delete_effect(effect_id: int) -> tuple[bool, dict[str, list]]:
        evidences = list(Evidence.select().where(Evidence.effect_id == effect_id))
        violations = {}
        if evidences:
            violations['Evidences'] = evidences
        return len(violations) == 0, violations

    @staticmethod
    def can_delete_subtype(subtype_id: int) -> tuple[bool, dict[str, list]]:
        evidences = list(Evidence.select().where(Evidence.subtype_id == subtype_id))
        violations = {}
        if evidences:
            violations['Evidences'] = evidences
        return len(violations) == 0, violations


class DatabaseOperations:
    """Database CRUD operations with transaction support."""

    @staticmethod
    def delete_segment(segment_name: str) -> None:
        """Delete a segment (must have no proteins or references)."""
        can_delete, violations = DeletionValidator.can_delete_segment(segment_name)
        if not can_delete:
            raise ForeignKeyViolationError('Segment', segment_name, violations)

        Segment.delete().where(Segment.name == segment_name).execute()
        Segment.clear_cache()

    @staticmethod
    def delete_protein(protein_id: int) -> None:
        """Delete a protein (must have no annotations or mutations)."""
        can_delete, violations = DeletionValidator.can_delete_protein(protein_id)
        if not can_delete:
            raise ForeignKeyViolationError('Protein', str(protein_id), violations)

        Protein.delete().where(Protein.id == protein_id).execute()

    @staticmethod
    def delete_reference(reference_id: int) -> None:
        """Delete a reference (must have no annotations or mappings)."""
        can_delete, violations = DeletionValidator.can_delete_reference(reference_id)
        if not can_delete:
            raise ForeignKeyViolationError('Reference', str(reference_id), violations)

        Reference.delete().where(Reference.id == reference_id).execute()
        Reference.clear_cache()

    @staticmethod
    def delete_annotation(annotation_id: int) -> None:
        """Delete an annotation (no dependencies)."""
        Annotation.delete().where(Annotation.id == annotation_id).execute()

    @staticmethod
    def delete_mutation(mutation_id: int) -> None:
        """Delete a mutation (must have no markers, will cascade delete mappings)."""
        can_delete, violations = DeletionValidator.can_delete_mutation(mutation_id)
        if not can_delete:
            raise ForeignKeyViolationError('Mutation', str(mutation_id), violations)

        # Delete mappings first (dependent records)
        Mapping.delete().where(Mapping.mutation_id == mutation_id).execute()
        # Delete mutation
        Mutation.delete().where(Mutation.id == mutation_id).execute()

    @staticmethod
    def delete_mapping(mapping_id: int) -> None:
        """Delete a mapping (no dependencies)."""
        Mapping.delete().where(Mapping.id == mapping_id).execute()

    @staticmethod
    def delete_marker(marker_id: int) -> None:
        """Delete a marker (must have no evidences)."""
        can_delete, violations = DeletionValidator.can_delete_marker(marker_id)
        if not can_delete:
            raise ForeignKeyViolationError('Marker', str(marker_id), violations)

        # Delete M2M associations
        Marker.mutations.get_through_model().delete().where(Marker.mutations.get_through_model().marker_id == marker_id).execute()
        # Delete marker
        Marker.delete().where(Marker.id == marker_id).execute()
        Marker.clear_cache()

    @staticmethod
    def delete_paper(paper_id: int) -> None:
        """Delete a paper (must have no evidences)."""
        can_delete, violations = DeletionValidator.can_delete_paper(paper_id)
        if not can_delete:
            raise ForeignKeyViolationError('Paper', str(paper_id), violations)

        Paper.delete().where(Paper.id == paper_id).execute()

    @staticmethod
    def delete_effect(effect_id: int) -> None:
        """Delete an effect (must have no evidences)."""
        can_delete, violations = DeletionValidator.can_delete_effect(effect_id)
        if not can_delete:
            raise ForeignKeyViolationError('Effect', str(effect_id), violations)

        Effect.delete().where(Effect.id == effect_id).execute()

    @staticmethod
    def delete_subtype(subtype_id: int) -> None:
        """Delete a subtype (must have no evidences)."""
        can_delete, violations = DeletionValidator.can_delete_subtype(subtype_id)
        if not can_delete:
            raise ForeignKeyViolationError('Subtype', str(subtype_id), violations)

        Subtype.delete().where(Subtype.id == subtype_id).execute()

    @staticmethod
    def delete_evidence(evidence_id: int) -> None:
        """Delete an evidence (no dependencies)."""
        Evidence.delete().where(Evidence.id == evidence_id).execute()

    @staticmethod
    def remove_marker_mutation_association(marker_id: int, mutation_id: int) -> None:
        """Remove a mutation from a marker (M2M)."""
        through_model = Marker.mutations.get_through_model()
        through_model.delete().where((through_model.marker_id == marker_id) & (through_model.mutation_id == mutation_id)).execute()
