from PySide6.QtCore import Qt
from PySide6.QtWidgets import QMenu, QTreeWidgetItem

from flumut_db_editor import database_operations
from flumut_db_editor.database_operations import ForeignKeyViolationError
from flumut_db_editor.gui.dialogs import (
    ConfirmationDialog,
    ErrorDialog,
    ForeignKeyViolationDialog,
    SuccessNotification,
)


class CrudMixin:
    """Context menu and delete workflow shared by table and tree tabs.

    Subclasses provide ``load_data()``, ``handle_edit()``, ``handle_delete()``
    and a ``get_selected_*`` helper.
    """

    def show_context_menu(self, pos):
        """Display the Edit / Delete / Refresh menu and dispatch the choice."""
        menu = QMenu(self)
        edit_action = menu.addAction('Edit')
        delete_action = menu.addAction('Delete')
        menu.addSeparator()
        refresh_action = menu.addAction('Refresh')

        action = menu.exec(self.mapToGlobal(pos))
        if action == edit_action:
            self.handle_edit()
        elif action == delete_action:
            self.handle_delete()
        elif action == refresh_action:
            self.on_refresh()

    def on_refresh(self):
        """Reload data from the database."""
        self.load_data()

    def handle_edit(self):
        """Edit the selected item - to be overridden by subclasses."""
        raise NotImplementedError

    def handle_delete(self):
        """Delete the selected item - to be overridden by subclasses."""
        raise NotImplementedError

    def confirm_and_delete(self, instance) -> bool:
        """Confirm with the user, delete ``instance``, and report the outcome.

        Returns True only if the record was deleted.
        """
        item_type = type(instance).__name__
        if not ConfirmationDialog.ask(
            self,
            f'Delete {item_type}?',
            f'Are you sure you want to delete this {item_type}?',
            'This action cannot be undone.',
        ):
            return False

        try:
            database_operations.delete(instance)
        except ForeignKeyViolationError as e:
            ForeignKeyViolationDialog.show_violation(self, e.item_type, e.item_name, e.violations)
            return False
        except Exception as e:
            ErrorDialog.show_error(self, 'Delete Failed', f'Failed to delete {item_type}.', str(e))
            return False

        SuccessNotification.show_success(self, f'{item_type} deleted successfully.')
        return True


class TableCrudMixin(CrudMixin):
    """CRUD mixin for table-based tabs (BaseTableTab)."""

    def get_selected_row(self) -> int | None:
        """Return the selected row index, or None (with a prompt) if none."""
        row = self.table.currentRow()
        if row < 0:
            ErrorDialog.show_error(self, 'No Selection', 'Please select an item first.')
            return None
        return row


class TreeCrudMixin(CrudMixin):
    """CRUD mixin for tree-based tabs (HierarchicalTab).

    Tree items carry their model instance in ``Qt.ItemDataRole.UserRole``.
    """

    def get_selected_item(self) -> QTreeWidgetItem | None:
        """Return the selected tree item, or None (with a prompt) if none."""
        item = self.tree.currentItem()
        if item is None:
            ErrorDialog.show_error(self, 'No Selection', 'Please select an item first.')
            return None
        return item

    def get_selected_instance(self):
        """Return the model instance stored on the selected tree item, if any."""
        item = self.get_selected_item()
        if item is None:
            return None
        return item.data(0, Qt.ItemDataRole.UserRole)
